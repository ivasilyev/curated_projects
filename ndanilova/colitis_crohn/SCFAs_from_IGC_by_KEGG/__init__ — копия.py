#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
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
import xlrd  # Required to parse Excel sheets
import requests
import bs4
import lxml  # Required for use with bs4
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
outputDir = "/data1/bio/projects/ndanilova/"
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
        return "\n".join(["\t".join([i] + self.as_dict[i]) for i in self.as_dict]) + "\n"
    @staticmethod
    def var_to_file(var_to_write, file_to_write):
        file = open(file_to_write, 'w')
        file.write(var_to_write)
        file.close()
    def dump(self, file):
        SampleDataArray.var_to_file(self.export(), file)


# Export sample data
wholeSampleDataArray = SampleDataArray.load_sampledata("/data2/bio/Metagenomes/SampleData/SAMPLEDATA_READS_SOLID_ILLUMINA_NO_HG19.txt")
processedSampleDataArray = SampleDataArray({i: wholeSampleDataArray.get(i) for i in sorted(list(set(wholeSampleDataArray.as_dict).intersection(groupDataDF["sample_name"].values.tolist())))})
sampleDataFileName = outputDir + "_".join(sorted(list(set(groupDataDF["group_id"].values.tolist())))) + ".sampledata"
processedSampleDataArray.dump(sampleDataFileName)


class KEGGCompoundSequencesRetriever:
    """
    This class processes KEGG web pages based on compound ID (CXXXXX)
    """
    def __init__(self, compound_id):
        self.c_id = compound_id
        self._headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/45.0.2454.85 Safari/537.36 OPR/32.0.1948.25'}
def parse_compound_page(self):
    """
    The function processes URL like 'http://www.genome.jp/dbget-bin/www_bget?cpd:<compound_id>'
    :return:
    Enzymes entries (X.X.X.X) list
    """
    _soup = bs4.BeautifulSoup(requests.get("http://www.genome.jp/dbget-bin/www_bget?cpd:" + c_id, headers=self._headers).content, "lxml")
    _


c_id = "C00246"
_headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/45.0.2454.85 Safari/537.36 OPR/32.0.1948.25'}
_soup = bs4.BeautifulSoup(requests.get("http://www.genome.jp/dbget-bin/www_bget?cpd:" + c_id, headers=_headers).content, "lxml")
table_soup = _soup.find_all("td", "fr2")[0]



row_soup = table_soup.find_all("tr")[13]
try:
    key_cell_soup = row_soup.find_all("th", "th20")[0]
    if key_cell_soup.text == "Enzyme":
        return [i.text for i in row_soup.find_all("td", "td20")[0].find_all("a")]
except IndexError:
    continue

key_cell.text == "Enzyme"

[i.text for i in value_cell.find_all("a")]

for row in _soup.find_all("tr"):
    key_cell = row.findall("th", "th20")
    value_cell = row.findall("td", "td20")



_soup.find_all("tr")
_soup.find_all("td", "td20")





processedGroupDataDF = baseGroupDataDF.loc[baseGroupDataDF["sample_name"].isin(list(sampleDataArray.as_dict))].append(pd.read_table("/data2/bio/Metagenomes/SampleData/GROUPDATA_SRR.tsv", sep='\t', header='infer', names=list(baseGroupDataDF), engine='python'))




list(set(sampleDataArray.as_dict).intersection(baseGroupDataDF["sample_name"].values.tolist()))
df.loc[df['channel'].isin(['sale','fullprice'])]


groupDataDF = groupDataDF.append(pd.read_table("/data2/bio/Metagenomes/SampleData/GROUPDATA_SRR.tsv", sep='\t', header='infer', names=list(groupDataDF), engine='python'))





file_to_list("/data2/bio/Metagenomes/SampleData/SAMPLEDATA_READS_SOLID_ILLUMINA_NO_HG19.txt")

samplesIDsDF["sample_name"].values.tolist() + controlsIDsDF["sample_name"].values.tolist()
filteredSampleDataList = [i for i in file_to_list("/data2/bio/Metagenomes/SampleData/SAMPLEDATA_READS_SOLID_ILLUMINA_NO_HG19.txt") if any(i.split("\t")[0] == j for j in groupDataDF["sample_name"].values.tolist() + groupDataDF["sample_name"].values.tolist())]

def var_to_file(var_to_write, file_to_write):
    file = open(file_to_write, 'w')
    file.write(var_to_write)
    file.close()

"""
# Look for Redis pod & service:
kubectl get pods --show-all

# Deploy if not present:
kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-pod.yaml && \
kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-service.yaml

# Copy master-pod.yaml & worker-job.yaml and deploy:
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ndanilova/master-pod.yaml

# Wait until master finish and create the job
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ndanilova/worker-job.yaml

# Look for pod
kubectl describe pod <NAME>

# Clean k8s
kubectl delete pod ndanilova-bwt-fp-queue && kubectl delete job ndanilova-bwt-fp-job

# Checkout
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker python3 /home/docker/scripts/verify_coverages.py -s /data1/bio/projects/ndanilova/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.sampledata -a /data/reference/IGC/index/igc_v2014.03_annotation.txt -m no_hg19_igc_v2014.03 -o /data1/bio/projects/ndanilova/Metagenomes/IGC
"""


requiredSampleDataFileName = outputDir + "_".join(sorted(list(set(samplesIDsDF["group_id"].values.tolist() + controlsIDsDF["group_id"].values.tolist())))) + ".sampledata"
var_to_file("\n".join(filteredSampleDataList), requiredSampleDataFileName)

# Now login to master node and




samplesIDsDF["sample_name"].values.tolist() + controlsIDsDF["sample_name"].values.tolist()
filteredSampleData = pd.read_table("/data2/bio/Metagenomes/SampleData/SAMPLEDATA_READS_SOLID_ILLUMINA_NO_HG19.txt", sep='\t', engine='python')





# Declare samples groups


groupsJson = {"colitis": {"remission": ["VZK5", "VZK8", "VZK28", "VZK30", "VZK86", "VZK92", "107VZK", "175VZK", "VZK 068", "VZK 76", "177VZK", "165VZK", "VZK55"],
                          "escalation": ["VZK94", "VZK7", "VZK9", "VZK12", "VZK13", "VZK14", "VZK15", "VZK19", "VZK20", "VZK21", "VZK24", "VZK29", "VZK32", "VZK34", "VZK45", "VZK47", "VZK51", "VZK52", "VZK54", "VZK56", "VZK78", "VZK93", "VZK81", "VZK83", "VZK84", "VZK85", "VZK87", "108VZK", "173VZK", "136VZK", "137VZK", "115VZK", "138VZK", "141VZK", "144VZK", "139VZK", "169VZK", "166VZK", "VZK 72", "147VZK", "164VZK", "167VZK", "161VZK", "148VZK", "156VZK", "99VZK", "146VZK", "VZK 59", "VZK 60", "VZK 63", "VZK 64", "VZK 66", "VZK 067", "VZK 69", "VZK 70", "VZK 71", "VZK 75", "VZK26"},
              "crohn": {"remission": []}





















def is_path_exists(path):
    try:
        os.makedirs(path)
    except OSError:
        pass


def make_cleanup(files_list):
    removed_files_list = []
    removed_directories_list = []
    for file_name in files_list:
        try:
            os.remove(file_name)
            removed_files_list.append(file_name)
        except FileNotFoundError:
            continue
        except IsADirectoryError:
            import shutil
            shutil.rmtree(file_name)
            removed_directories_list.append(file_name)
    if len(removed_files_list) > 0:
        print("Removed files: " + ', '.join(removed_file for removed_file in removed_files_list))
    if len(removed_directories_list) > 0:
        print("Removed directories: " + ', '.join(removed_dir for removed_dir in removed_directories_list))


def file_append(string, file_to_append):
    file = open(file_to_append, 'a+')
    file.write(string)
    file.close()


def external_route(input_direction_list, output_direction):
    process = subprocess.Popen(input_direction_list, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    if not output_direction:
        return re.sub('[\r\n]', '', output.decode("utf-8"))
    else:
        file_append(output.decode("utf-8"), output_direction)


def append_group_df(df, group_list, group_id):
    return df.append(pd.DataFrame({"sample_name": [re.sub('^0+', "", re.findall('(\d+)', re.sub('[ \r\n]', "", i))[0]) + re.findall('([A-Za-z]+)', re.sub('[ \r\n]', "", i))[0] for i in group_list], "group_id": [group_id,] * len(group_list)}), ignore_index=True)


def parse_base(xlsx):
    # Well, that was a something. Hope I'll get better data next time.
    col_names_list = ["sample_name", "group_id"]
    xl = pd.ExcelFile(xlsx)
    df = pd.DataFrame(columns=col_names_list)
    for sheet_id, col_name, group_id in zip(["ремиссия ЯК", "обострение ЯК", "обострение БК", "ремиссия БК"], ["# образца/Sample #", "# образца/Sample #", "# образца", "# образца"], ["colitis_rem", "colitis_esc", "crohn_esc", "crohn_rem"]):
        df1 = xl.parse(sheet_id)
        try:
            df = append_group_df(df, df1[col_name].values.tolist(), group_id)
        except KeyError:
            print([sheet_id, col_name, group_id], list(df1))
            raise
    return df.loc[:, col_names_list]


def make_groupdata():
    col_names_list = ["sample_name", "group_id"]
    # "colitis_rem", "colitis_esc", "crohn_esc", "crohn_rem" -> "vzk"
    dfs_dict = {"vzk": pd.DataFrame({"sample_name": baseDF["sample_name"].values.tolist(), "group_id": ["vzk",] * len(baseDF)}).loc[:, col_names_list]}
    # "colitis_rem", "colitis_esc" -> "colitis"; "crohn_esc", "crohn_rem" -> "crohn"
    df = baseDF.copy()
    df["group_id"] = df["group_id"].apply(lambda x: re.findall('(.*)_', x)[0])
    dfs_dict["disease"] = df.copy()
    del df
    # "colitis_rem", "crohn_rem" -> "rem"; "colitis_esc", "crohn_esc" -> "esc"
    df = baseDF.copy()
    df["group_id"] = df["group_id"].apply(lambda x: re.findall('_(.*)', x)[0])
    dfs_dict["escalation"] = df.copy()
    del df
    dfs_dict["disease_escalation"] = baseDF.copy()
    output_dict = {}
    for key in dfs_dict:
        output_file = outputDir + "srr_" + key + ".groupdata"
        dfs_dict[key].append(controlDF, ignore_index=True).to_csv(output_file, sep='\t', header=False, index=False)
        output_dict[key] = output_file
    return output_dict


def multi_core_queue(function_to_parallelize, queue):
    import multiprocessing
    pool = multiprocessing.Pool()
    output = pool.map(function_to_parallelize, queue)
    pool.close()
    pool.join()
    return output


def get_kegg_orthologs_df():
    # Wiki: Short-chain fatty acids (SCFAs), also referred to as volatile fatty acids (VFAs), are fatty acids with two to six carbon atoms: formate, acetate, propionate, butyrate, isobutyrate, valerate, isovalerate
    kegg_compounds_dict = {"formate": "C00058", "acetate": "C00033", "propionate": "C00163", "butyrate": "C00246", "isobutyrate": "C02632", "valerate": "C00803", "isovalerate": "C08262"}
    kegg_compounds_df = pd.DataFrame(data={"compound_name": list(kegg_compounds_dict), "compound_id": [kegg_compounds_dict[i] for i in kegg_compounds_dict]}).loc[:, ["compound_name", "compound_id"]]
    is_path_exists(outputDir + "kegg/")
    kegg_compound_ids_file = outputDir + "kegg/compound_ids.txt"
    kegg_orthology_ids_file_mask = outputDir + "kegg/orthology_ids"
    kegg_compounds_df["compound_id"].to_csv(kegg_compound_ids_file, sep='\t', header=False, index=False)
    make_cleanup([outputDir + "kegg/kegg_products2orthologs.log"])
    external_route(["python3", keggScriptDir + "kegg_products2orthologs.py", "-i", kegg_compound_ids_file, "-n", "-o", kegg_orthology_ids_file_mask], outputDir + "kegg/kegg_products2orthologs.log")
    return pd.read_table(kegg_orthology_ids_file_mask + "_detailed_dataset.tsv", sep='\t', header=0, engine="python")


def mp_grep_annotation4gene_name(keyword):
    # return list(set([i.split('\t')[1] for i in subprocess.getoutput("grep -F " + keyword + " " + referenceAnnotationFile).split('\n') if len(i) > 0]))
    return [i.split("\t") for i in subprocess.getoutput("grep -F " + keyword + " " + referenceAnnotationFile).split('\n')]


def var_to_file(var_to_write, file_to_write):
    file = open(file_to_write, 'w')
    file.write(var_to_write)
    file.close()



    # The sequence header does not contain required substrings itself, however, the annotation does. See http://meta.genomics.cn/meta/dataTools
reference_annotation_genes_mp_buffer = [j for i in multi_core_queue(mp_grep_annotation4gene_name, pd.unique(keggOrthologsDF["orthology_id"].values).tolist()) for j in i]
reference_annotation_genes_df = pd.DataFrame([i for i in reference_annotation_genes_mp_buffer if len(i[0]) > 0], columns=["Gene ID", "Gene Name", "Gene Length", "Gene Completeness Status", "Cohort Origin", "Taxonomic Annotation(Phylum Level)", "Taxonomic Annotation(Genus Level)", "KEGG Annotation", "eggNOG Annotation", "Sample Occurence Frequency", "Individual Occurence Frequency", "KEGG Functional Categories", "eggNOG Functional Categories", "Cohort Assembled"])
# Make SCFA-related genes list
reference_annotation_genes_file = outputDir + "kegg/reference/igc_scfa_genes.txt"
reference_annotation_genes_df["Gene Name"].to_csv(reference_annotation_genes_file, sep='\t', header=False, index=False)
# Select only sequences with headers relevant to SCFA
make_cleanup([outputDir + "kegg/reference/igc_scfa.fasta", outputDir + "kegg/reference/keywords2fasta.log", outputDir + "kegg/reference/index", outputDir + "kegg/reference/cook_the_reference.log"])
external_route(["python3", keggScriptDir + "keywords2fasta.py", "-f", referenceSequenceFile, "-k", reference_annotation_genes_file, "-i", "-o", outputDir + "kegg/reference/igc_scfa.fasta"], outputDir + "kegg/reference/keywords2fasta.log")
# Create indexes
external_route(["python3", nBeeScriptDir + "cook_the_reference.py", "-i", outputDir + "kegg/reference/igc_scfa.fasta", "-o", outputDir + "kegg/reference/index"], outputDir + "kegg/reference/cook_the_reference.log")
is_path_exists(outputDir + "kegg/mapping/")
make_cleanup(subprocess.getoutput("ls -t " + outputDir + "kegg/mapping/*.sampledata").split("\n"))
sample_data_file_name = outputDir + "kegg/mapping/igc_scfa_no_hg19.sampledata"
# Combine filtered sample data to process. It may contain 2 or 3 columns, so it's better to not use Pandas
var_to_file("\n".join([[j for j in subprocess.getoutput("grep -F " + i + " " + filteredSampleDataFile).split("\n") if j.split("\t")[0] == i][0] for i in baseDF.append(controlDF, ignore_index=True)["sample_name"].values.tolist() if len(i) > 0]), sample_data_file_name)
# Mapping queue with 9 nodes
nodes_list = ",".join(["node" + str(i) for i in reversed(range(1, 10))])
external_route(["python3", nBeeScriptDir + "refdata2schedule.py", "-i", sample_data_file_name, "-r", outputDir + "kegg/reference/index/igc_scfa.refdata", "-n", nodes_list, "-m", "no_hg19", "-t", "32", "-o", outputDir + "kegg/mapping/"], outputDir + "kegg/mapping/refdata2schedule.log")
# Check for broken pipes at the end just in case
external_route(["python3", nBeeScriptDir + "verify_coverages.py", "-s", sample_data_file_name, "-a", outputDir + "kegg/reference/index/igc_scfa_annotation.txt -m no_hg19_igc_scfa", "-o", outputDir + "kegg/mapping"], outputDir + "kegg/mapping/verify_coverages.log")
sample_data_files_list = subprocess.getoutput("ls -t " + outputDir + "kegg/mapping/*.sampledata").split("\n")
# If number of *.sampledata files > 0, there were broken pipes and data reprocessing is required
if len(sample_data_files_list) > 1:
    # Get the latest sample data and start remapping
    external_route(["python3", nBeeScriptDir + "refdata2schedule.py", "-i", sample_data_files_list[0], "-r", outputDir + "kegg/reference/index/igc_scfa.refdata", "-n", nodes_list, "-m", "no_hg19", "-t", "32", "-o", outputDir + "kegg/mapping/"], outputDir + "kegg/mapping/refdata2schedule.log")











headers_df = reference_annotation_genes_df.merge(keggOrthologsDF, left_on="Gene Name", right_on="orthology_id", how='outer').dropna(axis=0, subset=["sequence_header"])
headers_df.loc[:, ["sequence_header", "Cohort Origin", "Cohort Assembled", "Taxonomic Annotation(Phylum Level)", "Taxonomic Annotation(Genus Level)", "ortholog_id", "enzyme_name", "eggNOG Annotation", "KEGG Functional Categories", "eggNOG Functional Categories", "compound_kegg_id", "compound_name"]]


    reference_annotation_genes_df = reference_annotation_genes_df.dropna(axis=1, subset="Gene ID")

    reference_annotation_genes_list = sorted([i for i in list(set([k for j in reference_annotation_genes_mp_buffer for k in j])) if len(i) > 0])


    reference_annotation_genes_file = outputDir + "kegg/reference/igc_scfa_genes.txt"
    # Make SCFA-related genes list
    reference_annotation_genes_df["Gene Name"].to_csv(reference_annotation_genes_file, sep='\t', header=False, index=False)
    # Select only sequences with headers relevant to SCFA
    external_route(["python3", keggScriptDir + "keywords2fasta.py", "-f", referenceSequenceFile, "-k", reference_annotation_genes_file, "-i", "-o", outputDir + "kegg/reference/igc_scfa.fasta"], outputDir + "kegg/reference/keywords2fasta.log")
    # Create indexes
    external_route(["python3", nBeeScriptDir + "cook_the_reference.py", "-i", outputDir + "kegg/reference/igc_scfa.fasta", "-o", outputDir + "kegg/reference/index"], outputDir + "kegg/reference/cook_the_reference.log")
    is_path_exists(outputDir + "kegg/mapping/")
    make_cleanup(subprocess.getoutput("ls -t " + outputDir + "kegg/mapping/*.sampledata").split("\n"))
    sample_data_file_name = outputDir + "kegg/mapping/igc_scfa_no_hg19.sampledata"
    # Combine filtered sample data to process. It may contain 2 or 3 columns, so it's better to not use Pandas
    var_to_file("\n".join([[j for j in subprocess.getoutput("grep -F " + i + " " + filteredSampleDataFile).split("\n") if j.split("\t")[0] == i][0] for i in baseDF.append(controlDF, ignore_index=True)["sample_name"].values.tolist() if len(i) > 0]), sample_data_file_name)
    # Mapping queue with 1-3, 5-9 nodes
    external_route(["python3", nBeeScriptDir + "refdata2schedule.py", "-i", sample_data_file_name, "-r", outputDir + "kegg/reference/index/igc_scfa.refdata", "-n", ",".join(["node" + str(i) for i in reversed(range(1, 10)) if i != 4]), "-m", "no_hg19", "-t", "32", "-o", outputDir + "kegg/mapping/"], outputDir + "kegg/mapping/refdata2schedule.log")
    # Check for broken pipes at the end just in case
    external_route(["python3", nBeeScriptDir + "verify_coverages.py", "-s", sample_data_file_name, "-a", outputDir + "kegg/reference/index/igc_scfa_annotation.txt -m no_hg19_igc_scfa", "-o", outputDir + "kegg/mapping"], outputDir + "kegg/mapping/verify_coverages.log")
    sample_data_files_list = subprocess.getoutput("ls -t " + outputDir + "kegg/mapping/*.sampledata").split("\n")
    # If number of *.sampledata files > 0, there were broken pipes and data reprocessing is required
    if len(sample_data_files_list) > 1:
        # Get the latest sample data and start remapping
        external_route(["python3", nBeeScriptDir + "refdata2schedule.py", "-i", sample_data_files_list[0], "-r", outputDir + "kegg/reference/index/igc_scfa.refdata", "-n", ",".join(["node" + str(i) for i in reversed(range(1, 10)) if i != 4]), "-m", "no_hg19", "-t", "32", "-o", outputDir + "kegg/mapping/"], outputDir + "kegg/mapping/refdata2schedule.log")


def mp_grep_annotation4gene_name(keyword):
    return list(set([i.split('\t')[1] for i in subprocess.getoutput("grep -F " + keyword + " " + referenceAnnotationFile).split('\n') if len(i) > 0]))


def improve_kegg_annotation():
    is_path_exists(outputDir + "kegg/group_stats")
    reference_annotation_genes_file = outputDir + "kegg/reference/igc_scfa_genes.txt"

    kegg_orthologs_df = pd.read_table(kegg_orthology_ids_file, sep='\t', header=0, engine='python')


reference_annotation_df = pd.read_table(referenceAnnotationFile, sep='\t', header="infer", names=["Gene ID", "Gene Name", "Gene Length", "Gene Completeness Status", "Cohort Origin", "Taxonomic Annotation(Phylum Level)", "Taxonomic Annotation(Genus Level)", "KEGG Annotation", "eggNOG Annotation", "Sample Occurence Frequency", "Individual Occurence Frequency", "KEGG Functional Categories", "eggNOG Functional Categories", "Cohort Assembled"], engine='python')
reference_annotation_df.loc[reference_annotation_df["KEGG Annotation"].isin(kegg_orthologs_list)]["Gene Name"].fillna("").values.tolist()


reference_annotation_orthologs_list = sorted(list(set([k for j in [i.split(";") for i in reference_annotation_df["KEGG Annotation"].fillna("").values.tolist() if len(i) > 0] for k in j])), key = len)

kegg_orthologs_df = pd.read_table(outputDir + "kegg/orthology_ids.txt", sep='\t', header=0, engine='python')
kegg_orthologs_df["KEGG Compound"].values.tolist()



is_path_exists(outputDir + "kegg/reference")


if __name__ == '__main__':
    pass

groupDataScriptDir = "/home/biouser/scripts/ivasilyev/statistical_tools/"
keggScriptDir = "/home/biouser/scripts/ivasilyev/Misc/"
nBeeScriptDir = "/home/biouser/scripts/ivasilyev/bowtie-tools/"
outputDir = "/data1/bio/projects/Danilova_Natalya/"
suffixesFile = "/data2/bio/Metagenomes/SampleData/SUFFIXES_CARD_MvirDB_TADB_VFDB_viromes.tsv"
filterFile = "/data2/bio/Metagenomes/SampleData/FILTER_non-taxon_words.tsv"
referenceSequenceFile = "/data/reference/IGC/760MetaHit_139HMP_368PKU_511Bac.fa.90_95"
referenceAnnotationFile = "/data/reference/IGC/IGC.annotation_OF.summary"
filteredSampleDataFile = "/data2/bio/Metagenomes/SampleData/SAMPLEDATA_READS_SOLID_ILLUMINA_NO_HG19.txt"

# Generate group data
baseDF = parse_base(outputDir + "1 База на статистику самый точный вариант от 14.05.2017 (1).xlsx")
controlDF = pd.read_table("/data2/bio/Metagenomes/SampleData/GROUPDATA_SRR.tsv", sep='\t', header='infer', names=list(baseDF), engine='python')

if not os.path.isdir(outputDir + "kegg/"):
    keggOrthologsDF = get_kegg_orthologs_df()


    if not os.path.isdir(outputDir + "kegg/mapping/"):
        process_custom_reference()
    if not os.path.isdir(outputDir + "kegg/group_stats"):







groupDataFilesDict = make_groupdata()
for groupDataKey in groupDataFilesDict:
    external_route(["python3", groupDataScriptDir + "suffixes2statistics.py", "-s", suffixesFile, "-g", groupDataFilesDict[groupDataKey], "-f", filterFile, "-t", "10", "-c", "srr", "-m", "u-test", "-r", "fdr_bh", "-i", "reference_id", "-v", "id_mapped_reads_per_million_sample_total_reads", "-a", "0.05", "-o", outputDir + groupDataKey], outputDir + groupDataKey + "_suffixes2statistics.log")






















from numpy import nan
import multiprocessing
from collections import Counter
from matplotlib import pyplot as plt


def sanitize_bool(column_name, substring):
    def convert_bool(string, sub_string):
        string = re.sub('[\r\n  ]', '', str(string))
        if len(string) == 0 or string == 'nan':
            return nan
        return sub_string in string
    return outputDF.loc[:, column_name].apply(lambda x: convert_bool(x, substring))


def make_cleanup(files_list):
    removed_files_list = []
    removed_directories_list = []
    for file_name in files_list:
        try:
            os.remove(file_name)
            removed_files_list.append(file_name)
        except FileNotFoundError:
            continue
        except IsADirectoryError:
            import shutil
            shutil.rmtree(file_name)
            removed_directories_list.append(file_name)
    if len(removed_files_list) > 0:
        print("Removed files: " + ', '.join(removed_file for removed_file in removed_files_list))
    if len(removed_directories_list) > 0:
        print("Removed directories: " + ', '.join(removed_dir for removed_dir in removed_directories_list))


def file_append(string, file_to_append):
    file = open(file_to_append, 'a+')
    file.write(string)
    file.close()


def external_route(input_direction_list, output_direction):
    process = subprocess.Popen(input_direction_list, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    if not output_direction:
        return re.sub('[\r\n]', '', output.decode("utf-8"))
    else:
        file_append(output.decode("utf-8"), output_direction)


def external_route_wrapper(args):
    print(external_route(*args))


def multi_core_queue(function_to_parallelize, queue):
    pool = multiprocessing.Pool()
    pool.map(function_to_parallelize, queue)
    pool.close()
    pool.join()


def is_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        print(exception)


#########
inputFile = '/data1/bio/projects/Danilova_Natalya/база.xlsx'

outputDir = '/'.join(inputFile.split('/')[:-1]) + '/'
inputDF = pd.read_excel(inputFile)

# header:
# '# образца/Sample #VZK', 'diagn', 'remis/obostr', 'кальпротектина, норма до   50 мг/кг/<50 мкг/г норма;mcg/g normal to 50', 'город/city-1, село, деревня/village-2', 'woman-1,man-2', 'age', 'ДЛЯ ЯК степень поражения/colon involvement тотальное/total-1; левостороннее/left-sided -2, проктит/proctitis-3', 'ДЛЯ БК 1.-илеит; 2-колит: толстая кишка; 3-илеоколит, 4-баугиневая заслонка;5-2,4; 6 -3,4', '1. редко рецид-1 раз в год/rare relapses - once a year; 2.часто рец-более 2 раз в год/frequent relapses - more than twice a year', 'в/выявлено/initial diagnosis -1 , повторное/relapse-2', 'Тяжесть заболевания ДЛЯ ЯК /Diseases severity:легкое течение/mild(1), средтяж/moderate (2), тяж/severe (3), ремиссия/remission(4)', 'стадии для БК: воспалительная 1; стриктуриющая 2; фистулизирующая 3;4-2,3;5-1,3; 6-1,2', 'боль в животе/abdominal pain-1-есть', 'вздутие/bloating-2', 'слабость/fatigue-3', 'частота стула/stool frequency-4', 'частота стула/stool frequency:1. < 4   2.>4  3.>6', 'кровь и слизь в стуле presence of blood or mucus in stool -5-есть/present, 6-нет/absebt', 'внекишечные проявления/extraintestinal symptoms внекишечные проявления артралгия-1; афтозный стоматит -2; 3-вторичный спондилоартрит ассоц с ЯК, 2-х сакроилеит; 4-энтезопатия(энтезит первого реберно-грудинного сочленения слева, пяточные шпоры) ; 5-pyoderma gangrenosum', 'рост/height', 'вес/body weight', 'ИМТ/BMI', 'Курение/Smoking 1-yes, 2- no\n1-да, 2-нет.', 'Срок установления 1.(меньше 1 года) 2.(от 1 до 5 лет), 3.(от 5 до 10)4.от 10', 'Наследственность/Family predisposition: не отягощена/no -1,отягощена/yes -2', 'Стероидозавис/Steroids dependence.-1 нет/no 2-да/yes', 'Стероидорезист 1 нет;2-да/Steroids resistance 1/no, 2 - yes', 'ДЛЯ БК индекс Беста/Best index (for CD)', 'ДЛЯ БК индекс Беста:1.<150 неактив ст,2. 150-300 баллов легкая;3. 301-450 среднетяж;4. >450 тяжелая атака', 'ДЛЯ ЯК частичное Мейо', 'Мейо/Mayo index, сумма баллов', 'индекс Мейо/Mayo Score 0—2 — ремиссия/минимальная активность заболевания / remission/minimal disease activity-1 \n3—5 — легкая форма ЯК/mild UC-2\n6—10 — среднетяжелая форма ЯК/moderate UC-3\n11—12 — тяжелая форма ЯК/severe UC-4', 'ягл: 1-отр; 2 пол', 'анализ кала на скр кровь: 1.отр; 2.пол', 'инвалидность', 'Hb-гемоглобин, g/l', 'эритроциты Er x 10^12/L(3.70 - 5.20)', 'лейкоциты leu * 10 9', 'СОЭ/ESR(men 3 - 10mmmol per hour.Woman5 - 15 mmmol per hour)', 'п/stab(1- 6%)', 'с-я/(47-72%)segmented', 'эозинофилы/eosinophils (0,5-5%)', 'моноциты/monocytes (3-11%)', 'лимфоциты/lymphocytes(19-37%)', 'тром/platelets x 10^9/L', 'альбумин/albumin %(35-55)', 'общий белок г/л/ total protein g/l (64-83)', 'AST, U/L(5-34)', 'ALT, U/L(0-55)', 'Bilirubin, Total (3,4-20,5 umol/L)', 'Щф/alkaline phosphatase U/L(40-150)', 'GGT, U/L(9-64)', 'Glu(4,1-6,2) ммоль/л, mmol/l\nmmol/l', "ПТИ(70-130%)\nQuick's value", 'АСЛО/antistreptolysin<200  mE / ml', 'СРБ /C-reactive protein', 'ЦИК(N62-136)', 'IgA(1.1-3.5 мг/мл)', 'IgA(1.1-3.5 мг/мл).1', 'IgM(0.7-2.5 мг/мл)', 'IgG(6.65-16.45мг/мл)', 'РФ-ревматоидный фактор/<8-20\nrheumatoid factor', 'Al%', 'α1%', 'α2%', 'β%', 'γ %', '5 -АСК -1', 'гормоны принимал-1', 'азатиоприн 1-принимал', 'биологич терапия принимал-1'

outputDF = inputDF[inputDF.loc[:, '# образца/Sample #VZK'].notnull()].rename(columns={'# образца/Sample #VZK': 'sample_name',
                                                                                      'diagn': 'diagnosis',
                                                                                      'remis/obostr': 'is_in_remission',
                                                                                      'кальпротектина, норма до   50 мг/кг/<50 мкг/г норма;mcg/g normal to 50': 'calprotectin_mg_per_kg',
                                                                                      'город/city-1, село, деревня/village-2': 'is_from_city',
                                                                                      'woman-1,man-2': 'is_woman',
                                                                                      'ДЛЯ ЯК степень поражения/colon involvement тотальное/total-1; левостороннее/left-sided -2, проктит/proctitis-3': 'colitis_degree',
                                                                                      'ДЛЯ БК 1.-илеит; 2-колит: толстая кишка; 3-илеоколит, 4-баугиневая заслонка;5-2,4; 6 -3,4': 'crohn_degree',
                                                                                      '1. редко рецид-1 раз в год/rare relapses - once a year; 2.часто рец-более 2 раз в год/frequent relapses - more than twice a year': 'is_2_or_more_relapses_per_year',
                                                                                      'в/выявлено/initial diagnosis -1 , повторное/relapse-2': 'is_relapse',
                                                                                      'Тяжесть заболевания ДЛЯ ЯК /Diseases severity:легкое течение/mild(1), средтяж/moderate (2), тяж/severe (3), ремиссия/remission(4)': 'colitis_stage',
                                                                                      'стадии для БК: воспалительная 1; стриктуриющая 2; фистулизирующая 3;4-2,3;5-1,3; 6-1,2': 'chron_stage',
                                                                                      'боль в животе/abdominal pain-1-есть': 'is_abdominal_pain',
                                                                                      'вздутие/bloating-2': 'is_bloating',
                                                                                      'слабость/fatigue-3': 'is_fatiqued',
                                                                                      'частота стула/stool frequency-4': 'stool_per_day',
                                                                                      'частота стула/stool frequency:1. < 4   2.>4  3.>6': 'stool_frequency',
                                                                                      'кровь и слизь в стуле presence of blood or mucus in stool -5-есть/present, 6-нет/absebt': 'is_blood_or_mucus_in_stool',
                                                                                      'внекишечные проявления/extraintestinal symptoms внекишечные проявления артралгия-1; афтозный стоматит -2; 3-вторичный спондилоартрит ассоц с ЯК, 2-х сакроилеит; 4-энтезопатия(энтезит первого реберно-грудинного сочленения слева, пяточные шпоры) ; 5-pyoderma gangrenosum': 'extraintestinal_symptoms_stage',
                                                                                      'рост/height': 'height',
                                                                                      'вес/body weight': 'weight',
                                                                                      'ИМТ/BMI': 'body_mass_index',
                                                                                      'Курение/Smoking 1-yes, 2- no\n1-да, 2-нет.': 'is_smoking',
                                                                                      'Срок установления 1.(меньше 1 года) 2.(от 1 до 5 лет), 3.(от 5 до 10)4.от 10': 'diagnosis_term',
                                                                                      'Наследственность/Family predisposition: не отягощена/no -1,отягощена/yes -2': 'is_family_predisposition',
                                                                                      'Стероидозавис/Steroids dependence.-1 нет/no 2-да/yes': 'is_steroid_addiction',
                                                                                      'Стероидорезист 1 нет;2-да/Steroids resistance 1/no, 2 - yes': 'is_steroid_resistance',
                                                                                      'ДЛЯ БК индекс Беста/Best index (for CD)': 'crohn_best_index',
                                                                                      'ДЛЯ БК индекс Беста:1.<150 неактив ст,2. 150-300 баллов легкая;3. 301-450 среднетяж;4. >450 тяжелая атака': 'crohn_best_stage',
                                                                                      'ДЛЯ ЯК частичное Мейо': 'colitis_partial_mayo_score',
                                                                                      'Мейо/Mayo index, сумма баллов': 'colitis_total_mayo_score',
                                                                                      'индекс Мейо/Mayo Score 0—2 — ремиссия/минимальная активность заболевания / remission/minimal disease activity-1 \n3—5 — легкая форма ЯК/mild UC-2\n6—10 — среднетяжелая форма ЯК/moderate UC-3\n11—12 — тяжелая форма ЯК/severe UC-4': 'colitis_mayo_stage',
                                                                                      'ягл: 1-отр; 2 пол': 'is_helminthiase',
                                                                                      'анализ кала на скр кровь: 1.отр; 2.пол': 'is_fecal_occult_blood',
                                                                                      'инвалидность': 'disability_group',
                                                                                      'Hb-гемоглобин, g/l': 'hemoglobin_g_per_l',
                                                                                      'эритроциты Er x 10^12/L(3.70 - 5.20)': 'red_blood_cells_mln_per_ul',
                                                                                      'лейкоциты leu * 10 9': 'white_blood_cells_bln_per_l',
                                                                                      'СОЭ/ESR(men 3 - 10mmmol per hour.Woman5 - 15 mmmol per hour)': 'erythrocyte_sedimentation_rate',
                                                                                      'п/stab(1- 6%)': 'banded_neutrophils_percentage',
                                                                                      'с-я/(47-72%)segmented': 'segmented_neutrophils_percentage',
                                                                                      'эозинофилы/eosinophils (0,5-5%)': 'eosinophils_percentage',
                                                                                      'моноциты/monocytes (3-11%)': 'monocytes_percentage',
                                                                                      'лимфоциты/lymphocytes(19-37%)': 'lymphocytes_percentage',
                                                                                      'тром/platelets x 10^9/L': 'platelets_bln_per_l',
                                                                                      'альбумин/albumin %(35-55)': 'albumin_percentage',
                                                                                      'общий белок г/л/ total protein g/l (64-83)': 'total_protein_g_per_l',
                                                                                      'AST, U/L(5-34)': 'aspartate_aminotransferase_units_per_l',
                                                                                      'ALT, U/L(0-55)': 'alanine_aminotransferase_units_per_l',
                                                                                      'Bilirubin, Total (3,4-20,5 umol/L)': 'bilirubin_umol_per_l',
                                                                                      'Щф/alkaline phosphatase U/L(40-150)': 'alkaline_phosphatase_units_per_l',
                                                                                      'GGT, U/L(9-64)': 'gamma-glutamyltransferase_units_per_l',
                                                                                      'Glu(4,1-6,2) ммоль/л, mmol/l\nmmol/l': 'glutamic_acid_mmol_per_l',
                                                                                      "ПТИ(70-130%)\nQuick's value": 'prothrombin_ratio',
                                                                                      'АСЛО/antistreptolysin<200  mE / ml': 'antistreptolysin_o_units_per_l',
                                                                                      'СРБ /C-reactive protein': 'c-reactive_protein_mg_per_l',
                                                                                      'ЦИК(N62-136)': 'circulating_immune_complexes_units_per_l',
                                                                                      'IgA(1.1-3.5 мг/мл)': 'immunoglobulin_a_mg_per_ml',
                                                                                      'IgM(0.7-2.5 мг/мл)': 'immunoglobulin_m_mg_per_ml',
                                                                                      'IgG(6.65-16.45мг/мл)': 'immunoglobulin_g_mg_per_ml',
                                                                                      'РФ-ревматоидный фактор/<8-20\nrheumatoid factor': 'rheumatoid factor_units_per_ml',
                                                                                      'Al%': 'apolipoprotein_percentage',
                                                                                      'α1%': 'alpha1-globulin_percentage',
                                                                                      'α2%': 'alpha2-globulin_percentage',
                                                                                      'β%': 'beta-globulin_percentage',
                                                                                      'γ %': 'gamma-globulin_percentage',
                                                                                      '5 -АСК -1': 'is_senescent_cell_antigen',
                                                                                      'гормоны принимал-1': 'is_hormone_therapy',
                                                                                      'азатиоприн 1-принимал': 'is_azathioprine_therapy',
                                                                                      'биологич терапия принимал-1': 'is_probiotic_therapy'})

outputDF['diagnosis'] = outputDF.loc[:, 'diagnosis'].apply(lambda var: str(var).strip().replace('uc', 'colitis').replace('cd', 'crohn'))

for columnName, true_code in zip(['is_in_remission', 'is_from_city', 'is_woman', 'is_2_or_more_relapses_per_year', 'is_relapse', 'is_abdominal_pain', 'is_fatiqued', 'is_blood_or_mucus_in_stool', 'is_smoking', 'is_family_predisposition', 'is_steroid_addiction', 'is_steroid_resistance', 'is_helminthiase', 'is_fecal_occult_blood', 'is_senescent_cell_antigen', 'is_hormone_therapy', 'is_azathioprine_therapy', 'is_probiotic_therapy'], ['r', '1', '1', '2', '2', '1', '3', '5', '1', '2', '2', '2', '2', '2', '1', '1', '1', '1']):
    outputDF[columnName] = sanitize_bool(columnName, true_code)

# Export the whole raw data
outputDF.to_csv(outputDir + '/raw_group_data.tsv', sep='\t', index=False)

# Read external groupdata and combine only columns with sample names and group IDs
combinedGroupsDataFrame = outputDF.loc[:, ['sample_name', 'diagnosis']].append(pd.read_table('/data2/bio/Metagenomes/SampleData/GROUPDATA_SRR.tsv', sep='\t', header='infer', names=['sample_name', 'diagnosis'], engine='python'), ignore_index=True)
combinedGroupsDataFrame.rename(columns={'diagnosis': 'group_id'}, inplace=True)

# Make linkers using the formatted list: [[<prefix with path>, <suffix with extension>]]
# The IGC data is too big, selecting KEGG entries is required
# Prepare KEGG Compound IDs
is_path_exists(outputDir + "igc_kegg")
keggCompoundsDict = {"formate": "C00058", "acetate": "C00033", "propionate": "C00163", "butyrate": "C00246"}
keggCompoundsFile = outputDir + "igc_kegg/kegg_compounds_id_list.txt"
keggOrthologsFile = outputDir + "igc_kegg/kegg_orthologs_id_list.tsv"
make_cleanup([keggCompoundsFile, keggOrthologsFile])
external_route(["echo", '\n'.join([keggCompoundsDict[i] for i in keggCompoundsDict])], keggCompoundsFile)
# Get the related dataframe with KEGG KOs
external_route(["python3", "/home/biouser/scripts/kegg_products2orthologs.py", "-i", keggCompoundsFile, "-o", keggOrthologsFile], None)

keggOrthologsDF = pd.read_table(keggOrthologsFile, sep='\t', header=0, engine='python')
# Get KO IDs and flatten the resulting 2D arrays to lists within the dict
keggOrthologsRawDict = {i: [l for k in [j.split(' ') for j in keggOrthologsDF.loc[keggOrthologsDF['KEGG Compound'] == keggCompoundsDict[i]].loc[:, 'KEGG Orthologs'].values.tolist()] for l in k] for i in keggCompoundsDict}
# Remove overlapping KO IDs
keggOrthologsMixedList = [j for i in keggOrthologsRawDict.values() for j in i]
keggOrthologsUniqueIDsDict = {i: [j for j in keggOrthologsRawDict[i] if keggOrthologsMixedList.count(j) == 1] for i in keggOrthologsRawDict}
# Selecting from the whole coverage data only strings with required KEGG KOs


def grep4kegg(sample_name):
    sample_coverage_file = "/data2/bio/Metagenomes/IGC/Statistics/" + sample_name + "_no_hg19_igc_v1.0_coverage_combined.txt"
    output_coverage_file = outputDir + "igc_kegg/sample_coverages/" + keggCompoundsName + "/" + sample_name + "_" + keggCompoundsName + ".txt"
    make_cleanup([output_coverage_file])
    # Get the header
    external_route(["head", "-n", "1", sample_coverage_file], output_coverage_file)
    # Get the containment
    [external_route(["grep", "-F", i, sample_coverage_file], output_coverage_file) for i in keggOrthologsUniqueIDsDict[keggCompoundsName]]
    print("Completed selecting data based on KEGG Orthology IDs from '" + sample_coverage_file + "' to '" + output_coverage_file + "'")


# for keggCompoundsName in keggCompoundsDict:
#     # keggCompoundsName = "acetate"
#     is_path_exists(outputDir + "igc_kegg/sample_coverages/" + keggCompoundsName)
#     # Using the path of IGC extracted statistics files as prefix
#     try:
#         multi_core_queue(grep4kegg, combinedGroupsDataFrame.loc[:, "sample_name"].values.tolist())
#     except MemoryError:
#         print("Out of memory during the samples coverage data stripping! Use the single-thread scenario!")
#         [grep4kegg(i) for i in combinedGroupsDataFrame.loc[:, "sample_name"].values.tolist()]
#     print("The samples coverage data stripped succesfully!")

# The resulting coverage file name: <prefix><sample name><suffix>
linkerDict = {'CARD': ["/data2/bio/Metagenomes/ABR/CARD/Statistics/", "_no_hg19_card_v1.1.5_coverage.txt"],
              'MvirDB': ["/data2/bio/Metagenomes/Toxins/MvirDB/Statistics/", "_no_hg19_mvirdb_v2012.04.28_coverage.txt"],
              'TADB': ["/data2/bio/Metagenomes/Toxins/TADB/Statistics/", "_no_hg19_tadb_v1.1_coverage.txt"],
              'VFDB': ["/data2/bio/Metagenomes/Toxins/VFDB/Statistics/", "_no_hg19_vfdb_v2017.04.13_coverage.txt"]}

# Append the linker with samples stripped coverage data
# for keggCompoundsName in keggCompoundsDict:
#     linkerDict["IGC_KOs_for_" + keggCompoundsName] = [outputDir + "igc_kegg/sample_coverages/" + keggCompoundsName + "/", "_" + keggCompoundsName + ".txt"]

multiCoreCMDList = []
for screenedDBName in linkerDict:
    screenedDBPrefix, screenedDBSuffix = linkerDict[screenedDBName]
    screenedDBDataFrame = combinedGroupsDataFrame.copy()
    screenedDBDataFrame['sample_name'] = screenedDBDataFrame.loc[:, 'sample_name'].apply(lambda var1: screenedDBPrefix + var1 + screenedDBSuffix)
    screenedDBDataFrame.to_csv(outputDir + screenedDBName + '.groupdata', sep='\t', header=False, index=False)
    # Build schedule for the p-values calculating
    # multiCoreCMDList.append((["python3", "/home/biouser/scripts/groupdata2statistics.py", "-g", outputDir + screenedDBName + '.groupdata', "-i", "reference_id", "-v", "id_total_relative_abundance", "-o", outputDir + "pvals/" + screenedDBName], None))
    # 'id_total_relative_abundance' IGC values are too small for analysis
    multiCoreCMDList.append((["python3", "/home/biouser/scripts/groupdata2statistics.py", "-g", outputDir + screenedDBName + '.groupdata', "-i", "reference_id", "-v", "id_mapped_reads_per_million_sample_total_reads", "-o", outputDir + "pvals/" + screenedDBName], None))

# Calculate p-values using the groupdata
# try:
#     multi_core_queue(external_route_wrapper, multiCoreCMDList)
# except MemoryError:
#     print("Out of memory during the samples p-values calculating! Use the single-thread scenario!")
[external_route_wrapper(i) for i in multiCoreCMDList]

nonTaxonsWordsList = [i.lower() for i in ["common", "type", "pilus", "fimbriae", "hypothetical", "protein", "family", "encoded", "aida", "enterobactin", "secretion", "effector", "polysaccharide", "binding", "capsule", "hemin", "fibronectin", "sigma", "agglutinin", "extracellular", "hydrogen", "hemorrhagic", "invasion", "curli", "hemolysin", "flagella", "putative", "multi", "alpha", "trigger", "fimbriae", "uptake", "pili", "cytolisin", "siderophore", "metallophrotease", "iron", "hcp", "regulation", "sensor", "system", "fibers", "synthesis", "transport", "siga", "opacity", "outer", "two", "component", "biosynthesis", "yersiniabactin", "streptococcal", "glycosyltransferases", "ehaa", "purine", "pilus", "erm", "ere", "aada", "rrna", "mdt", "rec", "sub", "tet", "aac", "full", "translocase", "tetracycline", "resistance", "full", "membrane", "multidrug", "isoleucine", "plasmid", "beta", "lactamase", "vpar", "elen", "ecol", "ec", "shel", "amet", "npun", "dacr", "rml", "stg", "gro", "ps", "nle", "elg", "eha", "short", "arginine", "serine", "chromosome", "permeases", "antigen", "ferrous", "alginate", "fimbrial", "pyoverdine", "trehalose", "sortase", "recycling", "choline", "proteins", "intracellular", "multiplication", "adhesion", "copper", "exporter", "complex", "vibriobactin", "utilization", "periplasmic", "chaperone", "cytochrome", "muturation", "enhanced", "pyridine", "dithiocarboxylic", "amino", "acid", "phytotoxin", "syringomycin", "glucosyltransferase", "glucosyltransferases", "afimbrial", "adhesin", "fimbrin", "like", "mannose", "specific", "precursor", "accessory", "factor", "cation", "efflux", "predicted", "homolog", "dimethyladenosine", "transferase", "undecaprenyl", "diphosphatase", "uncharacterized", "bacitracin", "macrolide", "export", "import", "symport", "simport", "flagellar", "ring", "ferric", "reductase", "oxidoreductase", "phosphoenolpyruvate", "phosphotransferase", "carbamate", "kinase", "multiphosphoryl", "transfer", "acetate", "kinase", "phosphopantetheine", "adenylyltransferase", "motor", "preprotein", "subunit", "gamma", "epsilon", "alanine", "racemase", "peptidoglycan", "associated", "probable", "aminoglycoside", "bacteriocin", "phage", "uncharacterized", "transcriptional", "synthase", "synthetase", "stable", "inheritance", "tetracyline", "gene", "detergent", "resistant", "inner", "lipoprotein", "autoinducer", "production", "surface", "layer", "glycoside", "acetyl", "anaerobic", "nitric", "small", "big", "domain", "assembly", "acetyl", "methil", "ethyl", "propyl", "fusaric", "streptomycin", "phosphotransferase", "macrolide", "proteobactin", "salmochelin", "specific", "enzyme", "antibiotic", "multiple", "nucleotide", "nucleoside", "hook", "pilu", "kallikrein", "manganese", "transporter", "phosphopantetheine", "adenylyltransferase", "toxin", "anatoxin", "ribosomal" , "methyltransferase", "basal", "body", "general", "function", "streptothricin", "acetyltransferase", "phosphatidylcholine", "hydroxamate", "dependent", "neomycin", "kanamycin", "bifunctional", "polymyxin", "operon", "leader", "peptidase", "integration", "host"]]


def list2pairwise_tuples(flat_list):
    tuples = []
    for list_name1 in flat_list:
        for list_name2 in flat_list:
            if list_name1 != list_name2 and sorted([list_name1, list_name2]) not in [sorted(sublist) for sublist in tuples]:
                tuples.append([list_name1, list_name2])
    return tuples


def df2counters_dict(input_df, save_prefix):
    if len(input_df) == 0:
        print("Empty dataframe to count words!")
        return (None,) * 2
    reference_ids_list = input_df["reference_id"].values.tolist()
    filtered_ids_list = ["_".join([j for j in re.findall('[A-Za-z]{4,}', i) if j.lower() not in nonTaxonsWordsList]) for i in reference_ids_list]
    keyword_lists_dict = {'species': [j[0] for j in [re.findall('[A-Z]{1}[a-z]{4,}_[a-z]{4,}', i) for i in filtered_ids_list] if len(j) > 0],
                            'genera': [j[0] for j in [re.findall('[A-Z]{1}[a-z]{4,}', i) for i in filtered_ids_list] if len(j) > 0],
                            'common_words': [k for j in [re.findall('[a-z]{4,}', i.lower()) for i in reference_ids_list] for k in j]}
    output_dict = {}
    for keyword_lists_key in keyword_lists_dict:
        keyword_list = keyword_lists_dict[keyword_lists_key]
        if len(reference_ids_list) > len(keyword_list):
            keyword_lists_dict[keyword_lists_key] = keyword_list + ["unknown",] * (len(reference_ids_list) - len(keyword_list))
        keyword_df = pd.DataFrame.from_dict(Counter(keyword_lists_dict[keyword_lists_key]), orient='index').reset_index()
        keyword_df = keyword_df.rename(columns={list(keyword_df)[0]: 'keyword', list(keyword_df)[1]: 'occurrences'}).sort_values(by='occurrences', ascending=False)
        keyword_df.to_csv("_".join([save_prefix, "counter", keyword_lists_key]) + ".tsv", sep='\t', index=False)
        output_dict[keyword_lists_key] = keyword_df
        if len(keyword_df) > keyWordTopRowsNumber + 1:
            keyword_top_df = keyword_df.head(keyWordTopRowsNumber).append(pd.Series({'keyword': 'other', 'occurrences': keyword_df.tail(len(keyword_df) - keyWordTopRowsNumber)['occurrences'].sum()}), ignore_index=True)
            keyword_top_df.to_csv("_".join([save_prefix, "counter", keyword_lists_key, "top", str(keyWordTopRowsNumber)]) + ".tsv", sep='\t', index=False)
            output_dict[keyword_lists_key] = keyword_top_df
    return output_dict


def counter_df2plot(df, plot_title, save_prefix):
    df = df.set_index("keyword")
    fig = plt.figure(figsize=(12, 5))
    rotated_pairs_df = df.apply(lambda x: 100 * x / float(x.sum())).transpose()
    ax1 = rotated_pairs_df.plot(kind='bar', stacked=True, label='Percentage', width=.99)
    plt.xticks(rotation=0)
    ax1.set_xlabel("Comparison pairs")
    ax1.set_ylabel("Percentage", color='black')
    ax1.tick_params('y', colors='black')
    ax1.set_ylim([0, 100])
    labels = []
    for j in rotated_pairs_df.columns:
        for i in rotated_pairs_df.index:
            if rotated_pairs_df.loc[i][j] > 0:
                labels.append(j)
            else:
                labels.append("")
    patches = ax1.patches
    for label, rect in zip(labels, patches):
        width = rect.get_width()
        if width > 0:
            ax1.text(rect.get_x() + width / 2., rect.get_y() + rect.get_height() / 2., label, ha='center',
                     va='center', size='xx-small')
    for label, rect in zip(df.sum(), patches):
        ax1.text(rect.get_x() + rect.get_width() / 2., 101, str(int(label)), ha='center', va='bottom')
    ax1.legend(loc=1).set_visible(False)
    plt.title(plot_title, y=1.05)
    plt.savefig(save_prefix + ".png", dpi=600)
    plt.close(fig)


keyWordTopRowsNumber = 10
# Specifying the groups order
mainGroupsList = sorted(pd.unique(combinedGroupsDataFrame['group_id']).tolist())
mainGroupsList = [mainGroupsList[-1]] + mainGroupsList[:-1]
# len(list) is always more than index by 1
visualizingListsDict = {'free-for-all': list2pairwise_tuples(mainGroupsList)}


# Counting majors
for screenedDBName in linkerDict:
    screenedDBAdjustedDF = pd.read_table(subprocess.getoutput("ls -d " + outputDir + "pvals/" + screenedDBName + "/adjusted*u-test*").split('\n')[0], sep='\t', header=0)
    is_path_exists(outputDir + "pvals/" + screenedDBName + "/counters")
    # Select data to visualize
    visualizingCountersDFsDict = {}
    for screenedDBAdjustedDFBoolColName in [screenedDBAdjustedDFColName for screenedDBAdjustedDFColName in list(screenedDBAdjustedDF) if '_is_rejected_by_fdr_bh' in screenedDBAdjustedDFColName]:
        screenedDBAdjustedDFBoolColNamePair = re.findall('p-value_(.*)_is_rejected_by_fdr_bh', screenedDBAdjustedDFBoolColName)[0]
        screenedDBAdjustedDFBoolColNamePairDF = screenedDBAdjustedDF.loc[screenedDBAdjustedDF[screenedDBAdjustedDFBoolColName] == True]
        screenedDBAdjustedDFBoolColNamePairCountersDFsDict = df2counters_dict(screenedDBAdjustedDFBoolColNamePairDF, outputDir + "pvals/" + screenedDBName + "/counters/" + "_".join([screenedDBName, screenedDBAdjustedDFBoolColNamePair]))
        for screenedDBAdjustedDFBoolColNameGroup in screenedDBAdjustedDFBoolColNamePair.split('_vs_'):
            screenedDBAdjustedDFBoolColNameGroupIDsDF = screenedDBAdjustedDF.loc[screenedDBAdjustedDF[screenedDBAdjustedDFBoolColName] == True].loc[screenedDBAdjustedDF['prevalent_' + screenedDBAdjustedDFBoolColNamePair] == screenedDBAdjustedDFBoolColNameGroup]
            screenedDBAdjustedDFBoolColNameGroupIDsCountersDFsDict = df2counters_dict(screenedDBAdjustedDFBoolColNameGroupIDsDF, outputDir + "pvals/" + screenedDBName + "/counters/" + "_".join([screenedDBName, screenedDBAdjustedDFBoolColNamePair, "prevalent", screenedDBAdjustedDFBoolColNameGroup]))
        # Combining data in the single dataframe
        for screenedDBAdjustedDFBoolColNamePairCountersDFsDictKey in screenedDBAdjustedDFBoolColNamePairCountersDFsDict:
            keyWordDataFrame = screenedDBAdjustedDFBoolColNamePairCountersDFsDict[screenedDBAdjustedDFBoolColNamePairCountersDFsDictKey]
            try:
                visualizingCountersDFsDict[screenedDBAdjustedDFBoolColNamePairCountersDFsDictKey].update({screenedDBAdjustedDFBoolColNamePair: keyWordDataFrame})
            except KeyError:
                visualizingCountersDFsDict[screenedDBAdjustedDFBoolColNamePairCountersDFsDictKey] = {screenedDBAdjustedDFBoolColNamePair: keyWordDataFrame}
    visualizingSlidesDict = {}
    # visualizingListsDictKey - control type
    is_path_exists(outputDir + "pvals/" + screenedDBName + "/counters/plots")
    for visualizingListsDictKey in visualizingListsDict:
        visualizingSlidesDict[visualizingListsDictKey] = {}
        # visualizingGroupPairList - comparison list
        for visualizingGroupPairList in visualizingListsDict[visualizingListsDictKey]:
            visualizingGroupPairString = '_vs_'.join(visualizingGroupPairList)
            # visualizingCountersDFsDictKey - counter type
            for visualizingCountersDFsDictKey in visualizingCountersDFsDict:
                # Reassign data by counter type and control type
                if visualizingGroupPairString in list(visualizingCountersDFsDict[visualizingCountersDFsDictKey]):
                    visualizingSlidesDictToUpdate = {visualizingGroupPairString: visualizingCountersDFsDict[visualizingCountersDFsDictKey]['_vs_'.join(visualizingGroupPairList)]}
                elif '_vs_'.join(visualizingGroupPairList[::-1]) in list(visualizingCountersDFsDict[visualizingCountersDFsDictKey]):
                    visualizingSlidesDictToUpdate = {visualizingGroupPairString: visualizingCountersDFsDict[visualizingCountersDFsDictKey]['_vs_'.join(visualizingGroupPairList[::-1])]}
                else:
                    continue
                try:
                    visualizingSlidesDict[visualizingListsDictKey][visualizingCountersDFsDictKey].update(visualizingSlidesDictToUpdate)
                except KeyError:
                    visualizingSlidesDict[visualizingListsDictKey][visualizingCountersDFsDictKey] = visualizingSlidesDictToUpdate
    for visualizingListsDictKey in visualizingSlidesDict:
        for visualizingCountersDFsDictKey in visualizingSlidesDict[visualizingListsDictKey]:
            visualizingCounterDF = pd.DataFrame()
            for visualizingGroupPairString in visualizingSlidesDict[visualizingListsDictKey][visualizingCountersDFsDictKey]:
                if len(visualizingCounterDF) == 0:
                    visualizingCounterDF = visualizingSlidesDict[visualizingListsDictKey][visualizingCountersDFsDictKey][visualizingGroupPairString].rename(columns={"occurrences": visualizingGroupPairString})
                else:
                    visualizingCounterDF = pd.merge(visualizingCounterDF, visualizingSlidesDict[visualizingListsDictKey][visualizingCountersDFsDictKey][visualizingGroupPairString].rename(columns={"occurrences": visualizingGroupPairString}), on="keyword", how='outer')
                counter_df2plot(visualizingCounterDF, "_".join([screenedDBName, visualizingListsDictKey, visualizingCountersDFsDictKey]), outputDir + "pvals/" + screenedDBName + "/counters/plots/" + "_".join([screenedDBName, visualizingListsDictKey, visualizingCountersDFsDictKey]))
