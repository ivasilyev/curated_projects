#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git pull
python3
"""

import os
import pandas as pd
from meta.scripts.Utilities import Utilities
from inicolaeva.klebsiella_infants.ProjectDescriber import ProjectDescriber
import gzip
from Bio import SeqIO


# Get the raw reads files
raw_reads_files_dir = "/data1/bio/190405_M01969_0041_000000000-C6B66/Conversion_shotgun/Klebsiella"
raw_reads_files_list = [i for i in Utilities.scan_whole_dir(raw_reads_files_dir) if
                        os.path.normpath(os.path.dirname(i)) == raw_reads_files_dir and i.endswith("_001.fastq.gz")]
# Split them into the two groups
raw_reads_dict = {
    i: sorted([j for j in raw_reads_files_list if "_{}_".format(i) in os.path.splitext(os.path.basename(j))[0]]) for i
    in ("R1", "R2")}
# Combine the dict into the pandas.DataFrame object
raw_sampledata_df = pd.DataFrame.from_dict(raw_reads_dict)
# Are reads files corresponding to each other?
assert all((raw_sampledata_df["R1"].str.replace("_R1_", "_R2_") == raw_sampledata_df["R2"]).values.tolist() + (
            raw_sampledata_df["R2"].str.replace("_R2_", "_R1_") == raw_sampledata_df["R1"]).values.tolist())
# Get the sample names from reads file names
raw_sampledata_df["sample_name"] = raw_sampledata_df["R1"].map(
    lambda x: Utilities.safe_findall("(.+)_S[0-9]{2}_R[1|2]_001.fastq.gz", os.path.basename(x)))
# Export sampledata
project_describer = ProjectDescriber()
raw_sampledata_file = os.path.join(project_describer.ROOT_DIR, "sample_data", "raw_reads.sampledata")

Utilities.dump_tsv(df=raw_sampledata_df, table_file=raw_sampledata_file, col_names=["sample_name", "R1", "R2"])

print(raw_sampledata_file)  # /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/raw_reads.sampledata
# Create more detailed sampledata
raw_sampledata_df["reads_files"] = raw_sampledata_df.loc[:, ["R1", "R2"]].apply(lambda x: ";".join(x), axis=1)
raw_sampledata_df["taxon"] = "Klebsiella pneumoniae"
pipeline_sampledata_file = os.path.join(project_describer.ROOT_DIR, "sample_data", "raw_reads_pipeline.sampledata")

Utilities.dump_tsv(df=raw_sampledata_df, table_file=pipeline_sampledata_file,
                   col_names=["sample_name", "reads", "taxon"])

print(pipeline_sampledata_file)
# /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/raw_reads_pipeline.sampledata


def get_reads_stats_from_fq_gz(raw_reads_file):
    with gzip.open(raw_reads_file, "rt") as f:
        seq_records = list(SeqIO.parse(f, "fastq"))
        f.close()
    return {"sample_file": raw_reads_file, "sample_reads_number": len(seq_records),
            "sample_reads_bp": sum([len(i) for i in seq_records])}


reads_stats_list = Utilities.single_core_queue(get_reads_stats_from_fq_gz, raw_sampledata_df["R1"].values.tolist())
reads_stats_df = pd.DataFrame(reads_stats_list)
# Illumina's PE reads always have same counts of base pairs and total reads
raw_sampledata_df["sample_reads_number"] = reads_stats_df["sample_reads_number"] * 2
raw_sampledata_df["sample_reads_bp"] = reads_stats_df["sample_reads_bp"] * 2
# Count the expected coverage according to https://www.genome.jp/kegg-bin/show_organism?org=kpm
raw_sampledata_df["reference_genome_refseq"] = "NC_016845.1"
raw_sampledata_df["reference_genome_bp"] = 5333942
raw_sampledata_df["expected_coverage"] = raw_sampledata_df["sample_reads_bp"] / raw_sampledata_df["reference_genome_bp"]
raw_sampledata_df["expected_coverage"] = raw_sampledata_df["expected_coverage"].apply(
    lambda x: "{0:.1f}x".format(x))
reads_statistics_file = os.path.join(project_describer.ROOT_DIR, "sample_data", "reads_statistics.tsv")

Utilities.dump_tsv(raw_sampledata_df, reads_statistics_file,
                   col_names=["sample_name", "taxon", "sample_reads_number", "reference_genome_refseq",
                              "reference_genome_bp", "sample_reads_bp", "expected_coverage"])

print(reads_statistics_file)
# /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/reads_statistics.tsv
