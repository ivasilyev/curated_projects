#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup (on WORKER node):
docker pull python:latest && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it python:latest bash

docker pull ivasilyev/env_25_ecoli_genes:latest && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ivasilyev/env_25_ecoli_genes:latest bash
pip3 install pandas xlrd requests bs4 lxml jinja2 pyyaml
python3

"""

import os
import re
import pandas as pd
import xlrd  # Required by pandas to parse Excel sheets
import requests
import bs4
import lxml  # Required by bs4
import multiprocessing
import json
from collections import Counter
import subprocess
import yaml


def ends_with_slash(string):
    if string.endswith("/"):
        return string
    else:
        return str(string + "/")


def is_path_exists(path):
    try:
        os.makedirs(path)
    except OSError:
        pass


# Manage output location
outputDir = "/data1/bio/projects/ndanilova/colitis_crohn/"
is_path_exists(outputDir)
# Move the Excel DB table into outputDir and parse it


def parse_raw_database(xlsx):
    # Well, that was a something. Hope we'll get better data next time.
    col_names_list = ["sample_name", "group_id"]
    xl = pd.ExcelFile(xlsx)
    df = pd.DataFrame(columns=col_names_list)
    def append_group_dataframe(dataframe, group_list, group_id):
        return dataframe.append(pd.DataFrame({"sample_name": [re.sub('^0+', "", re.findall('(\d+)', re.sub('[ \r\n]', "", i))[0]) + re.findall('([A-Za-z]+)', re.sub('[ \r\n]', "", i))[0] for i in group_list], "group_id": [group_id,] * len(group_list)}), ignore_index=True)
    for sheet_id, col_name, group_id in zip(["ремиссия ЯК", "обострение ЯК", "обострение БК", "ремиссия БК"], ["# образца/Sample #", "# образца/Sample #", "# образца", "# образца"], ["colitis_rem", "colitis_esc", "crohn_esc", "crohn_rem"]):
        df1 = xl.parse(sheet_id)
        try:
            df = append_group_dataframe(df, df1[col_name].values.tolist(), group_id)
        except KeyError:
            print([sheet_id, col_name, group_id], list(df1))
            raise
    return df.loc[:, col_names_list]


# Parse base groups
groupDataDF = parse_raw_database(outputDir + "1 База на статистику самый точный вариант от 14.05.2017 (1).xlsx")
# And add the control group
groupDataDF = groupDataDF.append(pd.read_table("/data2/bio/Metagenomes/SampleData/GROUPDATA_SRR.tsv", sep='\t', header='infer', names=list(groupDataDF), engine='python'))

class SampleDataArray:
    def __init__(self, dictionary):
        self.as_dict = dict(dictionary)
    @staticmethod
    def load_sampledata(file):
        try:
            return SampleDataArray({j[0]: j[1:] for j in [i.split('\t') for i in SampleDataArray.file_to_list(file) if len(i) > 0]})
        except IndexError:
            raise ValueError("Cannot parse sample data: " + file)
    @staticmethod
    def file_to_list(file):
        file_buffer = open(file, 'rU')
        output_list = [j for j in [re.sub('[\r\n]', '', i) for i in file_buffer] if len(j) > 0]
        file_buffer.close()
        return output_list
    def get(self, key):
        return self.as_dict.get(key)
    def add(self, dictionary):
        self.as_dict.update(dictionary)
    def export(self):
        return "\n".join(["\t".join([i] + self.as_dict[i]) for i in self.as_dict])
    @staticmethod
    def var_to_file(var_to_write, file_to_write):
        file = open(file_to_write, 'w')
        file.write(var_to_write + "\n")
        file.close()
    def dump(self, file):
        self.var_to_file(self.export(), file)


# Export sample and group data
wholeSampleDataArray = SampleDataArray.load_sampledata("/data2/bio/Metagenomes/SampleData/SAMPLEDATA_READS_SOLID_ILLUMINA_NO_HG19.txt")
processedSampleDataArray = SampleDataArray({i: wholeSampleDataArray.get(i) for i in sorted(list(set(wholeSampleDataArray.as_dict).intersection(groupDataDF["sample_name"].values.tolist())))})
sampleDataFileSuffix = "_".join(sorted(list(set(groupDataDF["group_id"].values.tolist()))))
sampleDataFileName = outputDir + sampleDataFileSuffix + ".sampledata"
processedSampleDataArray.dump(sampleDataFileName)
groupDataDF.to_csv(outputDir + sampleDataFileSuffix + ".groupdata", sep='\t', header=True, index=False)


class FASTA:
    """
    This class is an attempt to apply NCBI standards to single FASTA.
    Consumes one header followed by sequence.
    """
    def __init__(self, single_fasta):
        self._body = re.sub("\n+", "\n", single_fasta.replace('\r', ''))
        try:
            self.header = re.findall(">(.+)", self._body)[0].strip()
            self.sequence = "\n".join(self.chunk_string(re.sub("[^A-Za-z]", "", self._body.replace(self.header, "")), 70)).upper()
        except IndexError:
            raise ValueError("Cannot parse the header!")
        # Nucleotide sequence has only AT(U)GC letters. However, it may be also protein FASTA.
    @staticmethod
    def chunk_string(string, length):
        return [string[0 + i:length + i] for i in range(0, len(string), length)]
    def __len__(self):
        return len(self.sequence)
    def to_dict(self):
        return {self.header: self.sequence}
    def to_str(self):
        return "\n".join([">" + self.header, self.sequence])


def remove_empty_values(input_list):
    try:
        return list(filter(lambda x: len(x) > 0, input_list))
    except TypeError:
        print(input_list)
        return []


def multi_core_queue(function_to_parallelize, queue):
    pool = multiprocessing.Pool()
    output = pool.map(function_to_parallelize, queue)
    pool.close()
    pool.join()
    return output


class KEGGCompoundSequencesRetriever:
    """
    This class processes KEGG web pages based on compound ID (CXXXXX)
    """
    def __init__(self, compound_id):
        self.c_id = compound_id
    @staticmethod
    def load_page(url):
        return bs4.BeautifulSoup(requests.get(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/45.0.2454.85 Safari/537.36 OPR/32.0.1948.25'}).content, "lxml")
    def parse_compound_page(self):
        """
        The function processes URL like 'http://www.genome.jp/dbget-bin/www_bget?cpd:<Compound ID>'
        :return:
        Enzymes entries (X.X.X.X) list
        """
        _soup = self.load_page("http://www.genome.jp/dbget-bin/www_bget?cpd:" + self.c_id)
        try:
            table_soup = _soup.find_all("td", "fr2")[0]
        except IndexError:
            return ""  # Page may be invalid
        for row_soup in table_soup.find_all("tr"):
            try:
                key_cell_soup = row_soup.find_all("th", "th20")[0]
                if key_cell_soup.text == "Enzyme":
                    enzyme_soup = row_soup.find_all("td", "td20")[0]
                    return [i.text for i in enzyme_soup.find_all("a")]
            except IndexError:
                continue
    def parse_enzyme_page(self, e_id):
        """
        The function processes URL like 'http://www.genome.jp/dbget-bin/www_bget?ec:<Enzyme ID>'
        :param e_id: str
        Enzyme ID corresponding the template 'X.X.X.X'.
        :return: list
        Orthologs list
        # Dictionary {'Orthology': [orthologs_list], 'Genes': [genes_list]}
        """
        _soup = self.load_page("http://www.genome.jp/dbget-bin/www_bget?ec:" + e_id)
        try:
            table_soup = _soup.find_all("td", "fr2")[0]
        except IndexError:
            return ""
        out_dict = {"Orthology": [], "Genes": []}
        for row_soup in table_soup.find_all("tr"):
            try:
                if row_soup.find_all("th", "th21")[0].text == "Orthology":
                    orthologs_soup = row_soup.find_all("td", "td21")[0]
                    out_dict["Orthology"].extend([i.text for i in orthologs_soup.find_all("a")])
            except IndexError:
                pass
            try:
                if row_soup.find_all("th", "th20")[0].text == "Genes":
                    for genes_soup in row_soup.find_all("td", "td20")[0].find_all("table")[0].find_all("tr"):
                        organism_id = genes_soup.find_all("nobr")[0].text.strip().lower()
                        gene_id = genes_soup.find_all("a")[0].text.strip()
                        out_dict["Genes"].append(organism_id + gene_id)
            except IndexError:
                pass
        return out_dict["Orthology"]
    def parse_ortholog_page(self, k_id):
        """
        The function processes URL like 'http://www.genome.jp/dbget-bin/www_bget?ko:<Orthology ID>'
        :param k_id: str
            Orthology ID corresponding the template 'KXXXXX'
        :return: list
        List [<Organism id>:<Gene ID>]
        """
        _soup = KEGGCompoundSequencesRetriever.load_page("http://www.genome.jp/dbget-bin/www_bget?ko:" + k_id)
        try:
            table_soup = _soup.find_all("td", "fr4")[0]
        except IndexError:
            return ""
        genes_list = []
        for row_soup in table_soup.find_all("tr"):
            try:
                if row_soup.find_all("th", "th40")[0].text == "Genes":
                    for gene_soup in row_soup.find_all("td", "td40")[0].find_all("table")[0].find_all("tr"):
                        organism_id = gene_soup.find_all("nobr")[0].text.strip().lower()
                        gene_ids_list = [i.text.strip() for i in gene_soup.find_all("a")]
                        genes_list.extend([organism_id + i for i in gene_ids_list])
            except IndexError:
                pass
        return genes_list
    def parse_gene_page(self, g_id):
        """
        The function processes URL like 'http://www.genome.jp/dbget-bin/www_bget?<Gene ID>'
        :param g_id: Gene ID corresponding the template '<Organism ID>:<Entry ID>'
        :return: Single FASTA object
        """
        _soup = self.load_page("http://www.genome.jp/dbget-bin/www_bget?" + g_id)
        try:
            table_soup = _soup.find_all("td", "fr1")[0]
        except IndexError:
            return ""
        out_dict = {"g_id": g_id}
        for row_soup in table_soup.find_all("tr"):
            try:
                if row_soup.find_all("th", "th11")[0].text.strip() == "Definition":
                    out_dict["definition"] = row_soup.find_all("td", "td11")[0].text
            except IndexError:
                pass
            try:
                if row_soup.find_all("th", "th10")[0].text.strip() == "KO":
                    out_dict["k_id"] = row_soup.find_all("td", "td10")[0].find_all("div")[0].find_all("div")[0].find_all("a")[0].text.strip()
                    out_dict["enzyme_name"] = row_soup.find_all("td", "td10")[0].find_all("div")[0].find_all("div")[1].text.strip()
            except IndexError:
                pass
            try:
                if row_soup.find_all("th", "th11")[0].text.strip() == "Organism":
                    out_dict["organism_name"] = " ".join([i.strip() for i in row_soup.find_all("td", "td11")[0].find_all("div")[0].text.split()[1:]])
            except IndexError:
                pass
            try:
                if row_soup.find_all("th", "th11")[0].text.strip() == "NT seq":
                    out_dict["sequence"] = "".join([i.strip().upper() for i in row_soup.find_all("td", "td11")[0].text.split('\n')[1:]])
            except IndexError:
                pass
        out_dict = {i: out_dict.get(i) if out_dict.get(i) else "" for i in ["g_id", "organism_name", "k_id", "enzyme_name", "definition", "sequence"]}
        return FASTA(">" + " ".join([out_dict.get(i).strip() for i in ["g_id", "organism_name", "k_id", "enzyme_name", "definition"]]) + "\n" + out_dict["sequence"])
    def compound2orthologs_list(self):
        e_ids_list = remove_empty_values(self.parse_compound_page())
        return sorted(remove_empty_values(list(set([k for j in multi_core_queue(self.parse_enzyme_page, e_ids_list) for k in j]))))
    def compound2genes_list(self):
        e_ids_list = remove_empty_values(self.parse_compound_page())
        k_ids_list = sorted(remove_empty_values(list(set([k for j in multi_core_queue(self.parse_enzyme_page, e_ids_list) for k in j]))))
        return sorted(remove_empty_values(list(set([k for j in multi_core_queue(self.parse_ortholog_page, k_ids_list) for k in j]))))
    def compound2fasta(self):
        e_ids_list = remove_empty_values(self.parse_compound_page())
        k_ids_list = sorted(remove_empty_values(list(set([k for j in multi_core_queue(self.parse_enzyme_page, e_ids_list) for k in j]))))
        g_ids_list = sorted(remove_empty_values(list(set([k for j in multi_core_queue(self.parse_ortholog_page, k_ids_list) for k in j]))))
        return remove_empty_values(multi_core_queue(self.parse_gene_page, g_ids_list))
    def genes_list2fasta(self, *args):
        # Just like the previous, but the 'g_ids_list' is given separately
        e_ids_list = remove_empty_values(self.parse_compound_page())
        k_ids_list = sorted(remove_empty_values(list(set([k for j in multi_core_queue(self.parse_enzyme_page, e_ids_list) for k in j]))))
        g_ids_list = sorted(remove_empty_values(list(set([k for j in multi_core_queue(self.parse_ortholog_page, k_ids_list) for k in j]).intersection(set(args)))))
        return remove_empty_values(multi_core_queue(self.parse_gene_page, g_ids_list))


# Create retrieving objects
scfaCompoundsIDsDict = {"formate": "C00058", "acetate": "C00033", "propionate": "C00163", "butyrate": "C00246", "isobutyrate": "C02632", "valerate": "C00803", "isovalerate": "C08262"}
scfaJSON = {}
foundGenesIDsList = []
for scfaCompoundName in scfaCompoundsIDsDict:
    scfaJSONDict = {"compound_name": scfaCompoundName, "compound_id": scfaCompoundsIDsDict[scfaCompoundName]}
    scfaJSONDict["retriever"] = KEGGCompoundSequencesRetriever(scfaJSONDict["compound_id"])
    scfaJSONDict.update({"orthologs_list": scfaJSONDict["retriever"].compound2orthologs_list(), "genes_list": scfaJSONDict["retriever"].compound2genes_list()})
    foundGenesIDsList.extend(scfaJSONDict["genes_list"])
    scfaJSON.update({scfaCompoundName: scfaJSONDict})

SampleDataArray.var_to_file(json.dumps({i: {j: scfaJSON[i][j] for j in scfaJSON[i] if j != "retriever"} for i in scfaJSON}), outputDir + "SCFAs_from_KEGG.json")

# Compile all FASTA objects
foundGenesIDsCounter = Counter(foundGenesIDsList)
notOverlappingGenesIDsList = [i for i in foundGenesIDsCounter if foundGenesIDsCounter[i] == 1]
notOverlappingGenesFASTAsList = [scfaJSON[i]["retriever"].genes_list2fasta(*notOverlappingGenesIDsList) for i in scfaJSON]

# Export FASTA
referenceDir = "/data/reference/custom/SCFAs_from_KEGG/"
is_path_exists(referenceDir)
SampleDataArray.var_to_file("\n".join([i.to_str() for i in remove_empty_values([k for j in notOverlappingGenesFASTAsList for k in j])]),
                            referenceDir + "SCFAs_from_KEGG.fasta")

"""
# Reference indexing
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker:latest \
python3 /home/docker/scripts/cook_the_reference.py \
-i /data/reference/custom/SCFAs_from_KEGG/SCFAs_from_KEGG.fasta \
-o /data/reference/custom/SCFAs_from_KEGG/index

"""

# Create config chart
chartsDir = outputDir + "SCFAs_from_KEGG/charts/"
is_path_exists(chartsDir)
cfgDict = {"QUEUE_NAME": "ndanilova-bwt-sk-queue",
           "MASTER_CONTAINER_NAME": "ndanilova-bwt-sk-master",
           "JOB_NAME": "ndanilova-bwt-sk-job",
           "ACTIVE_NODES_NUMBER": 7,
           "WORKER_CONTAINER_NAME": "ndanilova-bwt-sk-worker",
           "SAMPLEDATA": sampleDataFileName,
           "REFDATA": "/data/reference/custom/SCFAs_from_KEGG/index/SCFAs_from_KEGG.refdata",
           "OUTPUT_MASK": "no_hg19",
           "OUTPUT_DIR": "/data2/bio/Metagenomes/custom/SCFAs_from_KEGG"}
cfgFileName = chartsDir + "config.yaml"
with open(cfgFileName, 'w') as cfgFile:
    yaml.dump(cfgDict, cfgFile, default_flow_style=False, explicit_start=False)


def external_route(*args):
    process = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    return output.decode("utf-8")


# Dump script
genFileName = chartsDir + "generator.py"
external_route("curl", "-fsSL", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/generator.py", "-o", genFileName)

# Create charts from templates
external_route("python3", genFileName, "-c", cfgFileName, "-m", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/master.yaml", "-w", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/worker.yaml", "-o", chartsDir)

"""
# Copy charts into '/ndanilova/colitis_crohn/SCFAs_from_KEGG/' dir and push updates

# Pipeline launch
# Look for dependencies - Redis pod and service (from MASTER node)
kubectl get pods --show-all

# Deploy the MASTER chart to create queue
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ndanilova/colitis_crohn/SCFAs_from_KEGG/master.yaml

# Wait until master finish and deploy the WORKER chart to create the pipeline job
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ndanilova/colitis_crohn/SCFAs_from_KEGG/worker.yaml

# View active nodes
kubectl describe pod ndanilova-bwt-sk-job- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/custom/SCFAs_from_KEGG/Statistics/*coverage.txt | wc -l) OF $(cat /data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod ndanilova-bwt-sk-queue && \
kubectl delete job ndanilova-bwt-sk-job

# Checkout (from WORKER node)
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker python3 \
/home/docker/scripts/verify_coverages.py 
-s /data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.sampledata \
-r /data/reference/custom/SCFAs_from_KEGG/index/SCFAs_from_KEGG.refdata \
-m no_hg19_SCFAs_from_KEGG -d -o /data2/bio/Metagenomes/custom/SCFAs_from_KEGG

"""
