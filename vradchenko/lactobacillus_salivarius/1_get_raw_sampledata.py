#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git pull
LC_ALL=C python3
"""

import os
import pandas as pd
from meta.scripts.Utilities import Utilities
from vradchenko.lactobacillus_salivarius.ProjectDescriber import ProjectDescriber

# Get the raw reads files
raw_reads_files_dir = ProjectDescriber.RAW_DATA_DIR
raw_reads_files_list = [i for i in Utilities.scan_whole_dir(raw_reads_files_dir) if i.endswith("_001.fastq.gz")]

# Split them into the two groups
raw_reads_list = []
for raw_reads_files_pair in Utilities.get_most_similar_word_pairs(raw_reads_files_list):
    # Illumina file names have template '[sample]_[sequence]_[lane]_[strand]_[number].fastq.gz'
    # E.g: '336g_S1_L001_R1_001.fastq.gz'
    sample_name = Utilities.safe_findall("(.+)_S[0-9]+_L[0-9]+_R[0-9]+_[0-9]+", os.path.basename(raw_reads_files_pair[0]))
    raw_reads_dict = dict(sample_name=sample_name)
    STRANDS = ("R1", "R2")
    for raw_reads_file in raw_reads_files_pair:
        for reads_strand in STRANDS:
            if "_{}_".format(reads_strand) in os.path.splitext(os.path.basename(raw_reads_file))[0]:
                raw_reads_dict[reads_strand] = raw_reads_file
    if all([raw_reads_dict.get(STRANDS[0]).replace("_{}_".format(STRANDS[0]), "_{}_".format(STRANDS[-1])) == raw_reads_dict.get(STRANDS[-1])] +
           [raw_reads_dict.get(STRANDS[-1]).replace("_{}_".format(STRANDS[-1]), "_{}_".format(STRANDS[0])) == raw_reads_dict.get(STRANDS[0])]):
        raw_reads_list.append(raw_reads_dict)
    else:
        print("The read pair is invalid: '{}'".format(raw_reads_files_pair))

raw_sampledata_df = pd.DataFrame(raw_reads_list)
# Add suffices
raw_sampledata_df["sample_name"] = raw_sampledata_df["sample_name"] + "_" + raw_sampledata_df["R1"].apply(
    lambda x: Utilities.safe_findall("\w+", os.path.basename(os.path.dirname(x)).split("_")[-1]))
raw_sampledata_df["raw_reads"] = raw_sampledata_df["R1"] + ";" + raw_sampledata_df["R2"]


def define_species(_sample_name: str):
    _SPECIES = {"Bacillus subtilis BZR 336g": 336, "Bacillus subtilis BZR 517": 517,
                "Lactobacillus salivarius": 1, "Lactobacillus curvatus": 2, "Lactobacillus heilongjiangensis": 8}
    first_digits = int(Utilities.safe_findall("^\d+", _sample_name))
    for k in _SPECIES:
        if first_digits == _SPECIES.get(k):
            return k
    print("Cannot define species: '{}'".format(_sample_name))
    return "_"


raw_sampledata_df["host"] = raw_sampledata_df["sample_name"].apply(define_species)


Utilities.dump_tsv(raw_sampledata_df, ProjectDescriber.SAMPLE_DATA_FILE, col_names=["sample_name", "raw_reads"])
print(ProjectDescriber.SAMPLE_DATA_FILE)
# /data1/bio/projects/vradchenko/lactobacillus_salivarius/sample_data/raw.sampledata
