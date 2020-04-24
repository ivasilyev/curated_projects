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
from shutil import copy2
from meta.scripts.Utilities import Utilities
from vradchenko.lactobacillus_salivarius.ProjectDescriber import ProjectDescriber

# Get the raw reads files
raw_reads_files_dir = ProjectDescriber.RAW_DATA_DIR
raw_reads_files_list = [i for i in Utilities.scan_whole_dir(raw_reads_files_dir) if i.endswith("_001.fastq.gz")]

# Split them into the two groups
STRANDS = ("R1", "R2")
raw_reads_list = []
for raw_reads_files_pair in Utilities.get_most_similar_word_pairs(raw_reads_files_list):
    # Illumina file names have template '[sample]_[sequence]_[lane]_[strand]_[number].fastq.gz'
    # E.g: '336g_S1_L001_R1_001.fastq.gz'
    sample_name = Utilities.safe_findall("(.+)_S[0-9]+_L[0-9]+_R[0-9]+_[0-9]+", os.path.basename(raw_reads_files_pair[0]))
    raw_reads_dict = dict(sample_name=sample_name)
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
    first_digits = Utilities.safe_findall("^\d+", _sample_name)
    if len(first_digits) > 0:
        first_digits = int(first_digits)
        for k in _SPECIES:
            if first_digits == _SPECIES.get(k):
                return k
    print("Cannot define species: '{}'".format(_sample_name))
    return "_"


raw_sampledata_df["host"] = raw_sampledata_df["sample_name"].apply(define_species)

Utilities.dump_tsv(raw_sampledata_df, ProjectDescriber.SAMPLE_DATA_FILE, col_names=["sample_name", "raw_reads"])
print(ProjectDescriber.SAMPLE_DATA_FILE)
# /data1/bio/projects/vradchenko/lactobacillus_salivarius/sample_data/raw.sampledata

# Prepare Sequence Read Archive table
SRA_TEMPLATE_COL_NAMES = ["biosample_accession", "library_ID", "title", "library_strategy", "library_source", "library_selection",
                          "library_layout", "platform", "instrument_model", "design_description", "filetype", "filename", "filename2",
                          "filename3", "filename4", "assembly", "fasta_file", "bioproject_accession"]

sra_dir = os.path.join(ProjectDescriber.ROOT_DIR, "sra")
os.makedirs(os.path.join(sra_dir, "reads"), exist_ok=True)
sra_df = pd.DataFrame()
sra_df["sample_name"] = raw_sampledata_df["sample_name"]
for read_strand, filename_ in zip(STRANDS, ["filename", "filename2"]):
    sra_df["{}_source".format(read_strand)] = raw_sampledata_df[read_strand]
    sra_df[filename_] = sra_df["sample_name"] + "[{}].".format(read_strand) + sra_df["{}_source".format(read_strand)].apply(lambda x: ".".join(os.path.basename(x).split(".")[1:]))
    sra_df["{}_target".format(read_strand)] = sra_df[filename_].apply(lambda x: os.path.join(sra_dir, "reads", x))
    #
    _ = [copy2(*i) for i in sra_df.loc[:, ["{}_source".format(read_strand), "{}_target".format(read_strand)]].values]

for sra_regular_col_name, sra_regular_value in zip(
        ["library_strategy", "library_source", "library_selection", "library_layout", "platform", "instrument_model", "filetype"],
        ["WGS", "GENOMIC", "RANDOM", "paired", "ILLUMINA", "Illumina MiSeq", "fastq"]):
    sra_df[sra_regular_col_name] = sra_regular_value

sra_df["design_description"] = raw_sampledata_df["R1"].apply(lambda x: os.path.dirname(x).split("_")[-1])
sra_df["library_ID"] = sra_df["sample_name"]
sra_df.set_index("sample_name", inplace=True)

submission_report_df = Utilities.load_tsv("https://raw.githubusercontent.com/ivasilyev/curated_projects/master/vradchenko/lactobacillus_salivarius/data/tables/ncbi/submission_report.tsv").set_index("sample_name")
sra_df = pd.concat([sra_df, submission_report_df.loc[:, ["BioSample", "BioProject"]]], axis=1, sort=False)
sra_df.rename(columns={"BioSample": "biosample_accession", "BioProject": "bioproject_accession"}, inplace=True)

biosample_attributes_df = Utilities.load_tsv("https://raw.githubusercontent.com/ivasilyev/curated_projects/master/vradchenko/lactobacillus_salivarius/data/tables/ncbi/biosample_attributes_microbe.tsv").set_index("*sample_name")
sra_df = pd.concat([sra_df, biosample_attributes_df.loc[:, ["*organism", "isolation_source"]]], axis=1, sort=False)
sra_df["title"] = "DNA Miseq-PE-WGS of " + sra_df["*organism"] + ": " + sra_df["isolation_source"]

sra_df = sra_df.assign(**{i: "" for i in SRA_TEMPLATE_COL_NAMES if i not in sra_df.columns})

Utilities.dump_tsv(sra_df, os.path.join(sra_dir, "sra.tsv"), col_names=SRA_TEMPLATE_COL_NAMES)
