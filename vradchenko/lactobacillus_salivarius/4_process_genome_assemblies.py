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
from vradchenko.lactobacillus_salivarius.ProjectDescriber import ProjectDescriber

INDEX_COL_NAME = "sample_name"


def count_fasta_statistics(fasta_file: str, sample_name: str = None):
    from Bio import SeqIO
    with open(fasta_file, mode="r", encoding="utf-8") as f:
        seq_records = list(SeqIO.parse(f, "fasta"))
        f.close()
    out = dict(fasta_file=fasta_file, fasta_sequences_number=len(seq_records),
               fasta_total_bp=sum([len(i) for i in seq_records]))
    if sample_name:
        out["sample_name"] = sample_name
    return out


# Process assemblies
blasted_data_df = Utilities.load_tsv(os.path.join(ProjectDescriber.SAMPLE_DATA_DIR, "BLASTed.sampledata"))
assembly_files = blasted_data_df["assembly_file"].values.tolist()
assembly_stats_df = pd.DataFrame(Utilities.multi_core_queue(count_fasta_statistics, assembly_files))
assembly_stats_df.rename(columns={i: i.replace("fasta", "assembly") for i in assembly_stats_df.columns}, inplace=True)
blasted_data_df = pd.concat([blasted_data_df.set_index("assembly_file"), assembly_stats_df.set_index("assembly_file")],
                            axis=1, sort=False)
blasted_data_df.index.names = ["assembly_file"]
blasted_data_df = blasted_data_df.reset_index()

# Process raw reads
sample_data_df = Utilities.load_tsv(ProjectDescriber.SAMPLE_DATA_FILE)
sample_data_df["raw_file"] = sample_data_df["raw_reads"].apply(lambda x: x.split(";")[0].strip())
raw_read_files = sample_data_df["raw_file"].values.tolist()
raw_read_stats_df = pd.DataFrame(Utilities.multi_core_queue(Utilities.get_reads_stats_from_fq_gz, raw_read_files))
raw_read_stats_df.rename(columns={i: i.replace("sample", "raw") for i in raw_read_stats_df.columns}, inplace=True)
raw_read_stats_df["raw_reads_number"] *= 2
raw_read_stats_df["raw_reads_bp"] *= 2
sample_data_df = pd.concat([i.set_index("raw_file") for i in (sample_data_df, raw_read_stats_df)],
                           axis=1, sort=False)
sample_data_df.index.names = ["raw_file"]
sample_data_df = sample_data_df.reset_index()

combined_statistics_df = pd.concat([i.set_index(INDEX_COL_NAME) for i in (sample_data_df, blasted_data_df)],
                                   axis=1, sort=False)
combined_statistics_df.index.names = [INDEX_COL_NAME]
numeric_col_names = [i for i in combined_statistics_df.columns if any(i.endswith(j) for j in ("_bp", "_number"))]
combined_statistics_df.fillna(0, inplace=True)
combined_statistics_df = combined_statistics_df.astype({i: int for i in numeric_col_names})

combined_statistics_df["assembled_reads_percentage"] = combined_statistics_df["assembly_total_bp"] * 100 / combined_statistics_df["raw_reads_bp"]
combined_statistics_df["expected_assembly_coverage"] = combined_statistics_df["raw_reads_bp"] / combined_statistics_df["reference_bp"]
combined_statistics_df["real_assembly_coverage"] = combined_statistics_df["raw_reads_bp"] * (combined_statistics_df["assembled_reads_percentage"] / 100) / combined_statistics_df["reference_bp"]
for coverage_col_name in ("expected_assembly_coverage", "real_assembly_coverage"):
    combined_statistics_df[coverage_col_name] = combined_statistics_df[coverage_col_name].apply(lambda x: "{0:.2f}x".format(x))

combined_statistics_df = combined_statistics_df.reset_index()
combined_statistics_file = os.path.join(ProjectDescriber.SAMPLE_DATA_DIR, "combined_assembly_statistics.tsv")

Utilities.dump_tsv(combined_statistics_df.drop(columns=["raw_reads", ]), combined_statistics_file,
                   col_names=[i for i in combined_statistics_df.columns if "file" not in i])

print(combined_statistics_file)
# /data1/bio/projects/vradchenko/lactobacillus_salivarius/sample_data/combined_assembly_statistics.tsv
# Copy the data