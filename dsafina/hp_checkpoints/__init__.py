#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ivasilyev/bwt_filtering_pipeline_worker:latest python3

"""

import subprocess
import os
import re
import pandas as pd
import re


def is_path_exists(path):
    try:
        os.makedirs(path)
    except OSError:
        pass


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


def file_to_list(file):
    file_wrapper = open(file, 'rU')
    output_list = [j for j in [re.sub('[\r\n]', '', i) for i in file_wrapper] if len(j) > 0]
    file_wrapper.close()
    return output_list


def list_to_file(header, list_to_write, file_to_write):
    header += "\n".join(str(i) for i in list_to_write if i is not None)
    file = open(file_to_write, 'w')
    file.write(header)
    file.close()


if __name__ == '__main__':
    pass

outputDir = "/data1/bio/projects/dsafina/hp_checkpoints/"
is_path_exists(outputDir)
# Create combined groupdata
combinedGroupsDataFrame = pd.read_table("/data2/bio/Metagenomes/SampleData/GROUPDATA_HP_checkpoints_raw.tsv", sep='\t', header=0)
combinedGroupsDataFrame = combinedGroupsDataFrame.append(pd.read_table("/data2/bio/Metagenomes/SampleData/GROUPDATA_SRR.tsv", sep='\t', header='infer', names=list(combinedGroupsDataFrame), engine='python'), ignore_index=True)
combinedGroupsFileName = outputDir + "srr_hp_checkpoints.groupdata"
combinedGroupsDataFrame.to_csv(combinedGroupsFileName, sep='\t', header=False, index=False)
# Create combined sampledata
combinedSampleDataDict = {i.split('\t')[0]: i.split('\t')[1:] for i in file_to_list("/data2/bio/Metagenomes/SampleData/SAMPLEDATA_READS_SOLID_ILLUMINA_NO_HG19.txt") if len(i) > 0}
list_to_file("", ["\t".join([j] + combinedSampleDataDict[j]) for j in [i for i in combinedSampleDataDict if i in combinedGroupsDataFrame.loc[:, "sample_name"].values.tolist()]], outputDir + "srr_hp_checkpoints.sampledata")
# Process suffixes
groupDataScriptDir = "/home/biouser/scripts/ivasilyev/statistical_tools/"
suffixesFile = "/data2/bio/Metagenomes/SampleData/SUFFIXES_CARD_MvirDB_TADB_VFDB_viromes.tsv"
filterFile = "/data2/bio/Metagenomes/SampleData/FILTER_non-taxon_words.tsv"
external_route(["python3", groupDataScriptDir + "suffixes2statistics.py", "-s", suffixesFile, "-g", combinedGroupsFileName, "-f", filterFile, "-t", "10", "-c", "srr", "C", "-m", "u-test", "-r", "fdr_bh", "-i", "reference_id", "-v", "id_mapped_reads_per_million_sample_total_reads", "-a", "0.05", "-o", outputDir], outputDir + "suffixes2statistics.log")
