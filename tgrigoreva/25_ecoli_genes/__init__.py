#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${DOCKER_IMAGE_NAME} python3

"""

import re
import multiprocessing
import os
import requests
import bs4
import subprocess
import pandas as pd
import yaml

# Create reference sequences and metadata using gene name in initial HTML query
genesNamesList = "Stx2, EhxA, STb, EspA, EspB, EspC, Cnf, Cfa, Iha, pap, papA, papC, papE, papF, Tir, Etp, KpsM, KpsT, FliC, IbeA, Tsh, IucD, TraT, IutA, espP, katP, ompA, ompT, iroN, iss, fyuA, uidA, uspA, cdtB, cvaC, ibeA".split(", ")


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


def dict2pd_series(dictionary):
    output = pd.Series()
    for key in dictionary:
        output.at[key] = dictionary[key]
    return output


class NuccoreSequenceRetriever:
    """
    This class performs NCBI Gene DB search. Consumes organism name (space-delimited) and gene name.
    """
    def __init__(self, species, gene):
        self._species = species.strip()
        self._gene = gene.strip()
        self._headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.9; rv:45.0) Gecko/20100101 Firefox/45.0'}
        self._search_url = "https://www.ncbi.nlm.nih.gov/gene/?term=%22" + self._species.replace(' ', '+') + "%22%5Borgn%5D+AND+" + self._gene
        self._search_soup = bs4.BeautifulSoup(requests.get(self._search_url, headers=self._headers).content, "lxml")
        self._row_soups_list = self._search_soup.find_all("tr", "rprt")
    def get_soup(self):
        return self._search_soup
    @staticmethod
    def parse_table_row(row_soup):
        d = {"Name": row_soup.find_all("td", "gene-name-id")[0].find_all("a")[0].text,
             "Gene ID": "".join(re.findall("ID: (\d+)", row_soup.find_all("td", "gene-name-id")[0].find_all("span", "gene-id")[0].text)),
             "Description": row_soup.find_all("td")[1].text,
             "Sequence ID": "".join(re.findall("(NC_\d+\.\d*)", row_soup.find_all("td")[2].text)),
             "Sequence Location": "".join(re.findall("\((\d+\.\.\d+)", row_soup.find_all("td")[2].text)),
             "Aliases": row_soup.find_all("td")[3].text}
        return {k: d[k].strip() for k in d}
    def _get_fasta(self, organism_id, coordinates_list):
        if not isinstance(coordinates_list, list):
            raise ValueError("Coordinates must be in a list")
        coordinates_list = [str(i) for i in coordinates_list]
        _soup = bs4.BeautifulSoup(requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + organism_id + "&rettype=fasta&retmode=text&seq_start=" + coordinates_list[0] + "&seq_stop=" + coordinates_list[-1], headers=self._headers).content, "lxml")
        _fasta = FASTA(re.sub("\n+", "\n", _soup.find_all("p")[0].text))
        return _fasta
    def query2fasta(self):
        output_dict = {"FASTAs_list": [], "annotations_series_list": []}
        for _soup in self._row_soups_list:
            _row_dict = self.parse_table_row(_soup)
            _locations_list = _row_dict["Sequence Location"].split("..")
            # Filtering expression
            if (self._gene.lower() in _row_dict["Name"].lower() or self._gene.lower() in _row_dict["Description"].lower() or self._gene.lower() in _row_dict["Aliases"].lower()) and len(_locations_list) == 2:
                fasta = self._get_fasta(_row_dict["Sequence ID"], _locations_list)
                output_dict["FASTAs_list"].append(fasta)
                output_dict["annotations_series_list"].append(pd.Series(dict2pd_series(_row_dict), name=fasta.header))
        return output_dict


# Get sequences and metadata
def multi_core_queue(function_to_parallelize, queue):
    pool = multiprocessing.Pool()
    output = pool.map(function_to_parallelize, queue)
    pool.close()
    pool.join()
    return output


def mp_gene_search(name):
    obj = NuccoreSequenceRetriever("Escherichia coli", name)
    return obj.query2fasta()


foundGenesDictsList = multi_core_queue(mp_gene_search, genesNamesList)


# Merge dictionaries
def gene_search_wrapper():
    output_dict = {i: [] for i in foundGenesDictsList[0]}
    for input_dict in foundGenesDictsList:
        if len(input_dict[list(input_dict)[0]]) > 0:
            for key in input_dict:
                output_dict[key].extend(input_dict[key])
    df = pd.concat(output_dict["annotations_series_list"], axis=1).transpose()
    df.index.names = ["former_id"]
    return {"FASTAs_list": [i.to_str() for i in output_dict["FASTAs_list"]], "annotations_series_list": df}


foundGenesDict = gene_search_wrapper()


# Dump reference
def list_to_file(header, list_to_write, file_to_write):
    header += "\n".join(str(i) for i in list_to_write if i is not None) + "\n"
    file = open(file_to_write, 'w')
    file.write(header)
    file.close()


referenceDir = "/data/reference/custom/25_ecoli_genes/"
os.makedirs(referenceDir, exist_ok=True)
list_to_file("", foundGenesDict["FASTAs_list"], referenceDir + "25_ecoli_genes.fasta")

"""
# Reference indexing
rm -rf /data/reference/custom/25_ecoli_genes/index
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker:latest \
python3 /home/docker/scripts/cook_the_reference.py \
-i /data/reference/custom/25_ecoli_genes/25_ecoli_genes.fasta \
-o /data/reference/custom/25_ecoli_genes/index

"""


# Add query metadata to annotation
def update_annotation():
    annotation_file_name = referenceDir + "index/25_ecoli_genes_annotation.txt"
    df = pd.read_table(annotation_file_name, sep='\t', header=0, engine='python').set_index("former_id")
    df = pd.concat([df, foundGenesDict["annotations_series_list"]], axis=1).reset_index().set_index("reference_id").sort_index(ascending=True)
    df.to_csv(annotation_file_name, sep='\t', index=True, header=True)


update_annotation()


# Create charts
def external_route(*args):
    process = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    return output.decode("utf-8")


outputDir = "/data1/bio/projects/tgrigoreva/25_ecoli_genes/"
chartsDir = outputDir + "charts/"
subprocess.getoutput("rm -rf " + chartsDir)
os.makedirs(chartsDir, exist_ok=True)
cfgDict = {"QUEUE_NAME": "tgrigoreva-bwt-25-queue",
           "MASTER_CONTAINER_NAME": "tgrigoreva-bwt-25-master",
           "JOB_NAME": "tgrigoreva-bwt-25-job",
           "ACTIVE_NODES_NUMBER": 9,
           "WORKER_CONTAINER_NAME": "tgrigoreva-bwt-25-worker",
           "SAMPLEDATA": "/data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.sampledata",
           "REFDATA": "/data/reference/custom/25_ecoli_genes/index/25_ecoli_genes.refdata",
           "OUTPUT_MASK": "no_hg19",
           "OUTPUT_DIR": "/data2/bio/Metagenomes/custom/25_ecoli_genes"}

# Dump config
cfgFileName = chartsDir + "config.yaml"
with open(cfgFileName, 'w') as file:
    yaml.dump(cfgDict, file, default_flow_style=False, explicit_start=False)

# Dump script
genFileName = chartsDir + "generator.py"
external_route("curl", "-fsSL", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/generator.py", "-o", genFileName)

# Create charts from templates
external_route("python3", genFileName, "-c", cfgFileName, "-m", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/master.yaml", "-w", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/worker.yaml", "-o", chartsDir)

"""
# Copy charts into '/master/bwt_filtering_pipeline/' dir and push updates

# Pipeline launch
# Look for Redis pod & service:
kubectl get pods --show-all

# Deploy if not present:
kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-pod.yaml && \
kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-service.yaml

# Deploy the MASTER chart to create queue
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/tgrigoreva/25_ecoli_genes/master.yaml

# Wait until master finish and deploy the WORKER chart to create the pipeline job
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/tgrigoreva/25_ecoli_genes/worker.yaml

# View active nodes
kubectl describe pod tgrigoreva-bwt-25-job- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/custom/25_ecoli_genes/Statistics/*coverage.txt | wc -l) OF $(cat /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod tgrigoreva-bwt-25-queue
kubectl delete job tgrigoreva-bwt-25-job

# Checkout (from WORKER node)
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker python3 \
/home/docker/scripts/verify_coverages.py -s /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints.sampledata \
-r /data/reference/custom/25_ecoli_genes/index/25_ecoli_genes.refdata \
-m no_hg19_25_ecoli_genes -d -o /data2/bio/Metagenomes/custom/25_ecoli_genes

"""
