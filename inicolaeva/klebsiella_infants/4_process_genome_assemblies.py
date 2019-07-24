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
from Bio import SeqIO

assembler_result_dir = "/data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline/04_spades"
assembly_files = [i for i in Utilities.scan_whole_dir(assembler_result_dir) if os.path.basename(i) == "contigs.fasta"]
assemblies_target_dir = "/data1/bio/projects/inicolaeva/klebsiella_infants/assemblies"

assemblies_dict = dict()
for assembly_file in assembly_files:
    sample_name = Utilities.safe_findall("(Kleb[0-9]+)", assembly_file)
    if "plasmid" in assembly_file:
        assembly_type = "plasmid"
    else:
        assembly_type = "genome"
    if assemblies_dict.get(sample_name) is None:
        assemblies_dict[sample_name] = dict()
    assemblies_dict[sample_name]["sample_name"] = sample_name
    assemblies_dict[sample_name]["{}_file".format(assembly_type)] = assembly_file
    seq_records_raw = list(SeqIO.parse(assembly_file, "fasta"))
    assemblies_dict[sample_name]["{}_assembly_contigs_number_raw".format(assembly_type)] = len(seq_records_raw)
    assemblies_dict[sample_name]["{}_assembly_bp_raw".format(assembly_type)] = sum([len(i) for i in seq_records_raw])
    # NCBI does not allow to submit sequences shorter than 200 nucleotides
    seq_records_valid = [i for i in list(SeqIO.parse(assembly_file, "fasta")) if len(i) > 200]
    assemblies_dict[sample_name]["{}_assembly_contigs_number_valid".format(assembly_type)] = len(seq_records_valid)
    assemblies_dict[sample_name]["{}_assembly_bp_valid".format(assembly_type)] = sum(
        [len(i) for i in seq_records_valid])
    SeqIO.write(seq_records_valid, os.path.join(assemblies_target_dir, "{}_{}.fna".format(sample_name, assembly_type)),
                "fasta")


INDEX_COL_NAME = "sample_name"
assemblies_statistics_df = pd.DataFrame(assemblies_dict.values()).set_index(INDEX_COL_NAME)
reads_statistics_file = "/data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/reads_statistics.tsv"
reads_statistics_df = Utilities.load_tsv(reads_statistics_file).set_index(INDEX_COL_NAME)
combined_statistics_df = pd.concat([reads_statistics_df, assemblies_statistics_df], axis=1, sort=False)
combined_statistics_df.index.names = [INDEX_COL_NAME]
numeric_col_names = [i for i in combined_statistics_df.columns if any(
    j in i for j in ("_assembly_contigs_", "_assembly_bp_"))]
combined_statistics_df.fillna(0, inplace=True)
combined_statistics_df = combined_statistics_df.astype({i: int for i in numeric_col_names})
# From NCBI template ('Template_GenomeBatch.11700383121d.xlsx'):
# The estimated base coverage across the genome, eg 12x.
# This can be calculated by dividing the number of bases sequenced by the expected genome size
# and multiplying that by the percentage of bases that were placed in the final assembly.
# More simply it is the number of bases sequenced divided by the expected genome size.
combined_statistics_df["genome_assembled_reads_percentage"] = combined_statistics_df["genome_assembly_bp_valid"] * 100 / combined_statistics_df["sample_reads_bp"]
combined_statistics_df["genome_assembly_coverage"] = combined_statistics_df["sample_reads_bp"] * (combined_statistics_df["genome_assembled_reads_percentage"] / 100) / combined_statistics_df["reference_genome_bp"]
combined_statistics_df["genome_assembly_coverage"] = combined_statistics_df["genome_assembly_coverage"].apply(lambda x: "{0:.2f}x".format(x))
combined_statistics_file = os.path.join(ProjectDescriber.ROOT_DIR, "sample_data", "combined_assembly_statistics.tsv")
combined_statistics_df.reset_index(inplace=True)

Utilities.dump_tsv(combined_statistics_df, combined_statistics_file,
                   col_names=[i for i in combined_statistics_df.columns if not i.endswith("file")])

print(combined_statistics_file)
# /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/combined_assembly_statistics.tsv
# Copied into the ./datasets directory
