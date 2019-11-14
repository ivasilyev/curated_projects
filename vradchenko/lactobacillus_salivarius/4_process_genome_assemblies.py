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
import re
import pandas as pd
from shutil import copy2
from datetime import datetime
from meta.scripts.Utilities import Utilities
from meta.scripts.ncbi_contamination_remover import ContaminationRemover
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
blasted_data_df["organism"] = blasted_data_df["strain"].apply(lambda x: " ".join(x.split(" ")[:2]))

blasted_data_df.rename(columns={i: "reference_{}".format(i) for i in blasted_data_df.columns if
                                all(j not in i for j in ["assembly", "reference", "sample"])}, inplace=True)
assembly_files = blasted_data_df["assembly_file"].values.tolist()

assembly_stats_df = pd.DataFrame(Utilities.multi_core_queue(Utilities.count_assembly_statistics, assembly_files))
assembly_stats_df.rename(columns={i: "assembly_{}".format(i) for i in assembly_stats_df.columns}, inplace=True)
blasted_data_df = pd.concat([blasted_data_df.set_index("assembly_file"), assembly_stats_df.set_index("assembly_file")],
                            axis=1, sort=False)
blasted_data_df.index.names = ["assembly_file"]
blasted_data_df.reset_index(inplace=True)

# Process raw reads
sample_data_df = Utilities.load_tsv(ProjectDescriber.SAMPLE_DATA_FILE)
sample_data_df["raw_file"] = sample_data_df["raw_reads"].apply(lambda x: x.split(";")[0].strip())
raw_read_files = sample_data_df["raw_file"].values.tolist()


def mp_count_raw_reads_statistics(reads_file):
    return Utilities.count_raw_reads_statistics(reads_file=reads_file, type_="fastq_gz")


raw_read_stats_df = pd.DataFrame(Utilities.multi_core_queue(mp_count_raw_reads_statistics, raw_read_files))

raw_read_stats_df.rename(columns={i: "raw_{}".format(i) for i in raw_read_stats_df.columns}, inplace=True)
raw_read_stats_df["raw_reads_number"] *= 2
raw_read_stats_df["raw_total_reads_bp"] *= 2
sample_data_df = pd.concat([i.set_index("raw_file") for i in (sample_data_df, raw_read_stats_df)],
                           axis=1, sort=False)
sample_data_df.index.names = ["raw_file"]
sample_data_df.reset_index(inplace=True)

combined_statistics_df = pd.concat([i.set_index(INDEX_COL_NAME) for i in (sample_data_df, blasted_data_df)],
                                   axis=1, sort=False)
combined_statistics_df.index.names = [INDEX_COL_NAME]
numeric_col_names = [i for i in combined_statistics_df.columns if any(i.endswith(j) for j in ("_bp", "_number"))]
combined_statistics_df.fillna(0, inplace=True)
combined_statistics_df = combined_statistics_df.astype({i: int for i in numeric_col_names})

combined_statistics_df["assembled_reads_percentage"] = combined_statistics_df["assembly_total_contigs_bp"] * 100 / combined_statistics_df["raw_total_reads_bp"]
combined_statistics_df["expected_assembly_coverage"] = combined_statistics_df["raw_total_reads_bp"] / combined_statistics_df["reference_bp"]
combined_statistics_df["real_assembly_coverage"] = combined_statistics_df["raw_total_reads_bp"] * (combined_statistics_df["assembled_reads_percentage"] / 100) / combined_statistics_df["reference_bp"]
for coverage_col_name in ("expected_assembly_coverage", "real_assembly_coverage"):
    combined_statistics_df[coverage_col_name] = combined_statistics_df[coverage_col_name].apply(lambda x: "{0:.2f}x".format(x))


def define_strain_name(s: str):
    _OWNER_PREFIX = "VRA"
    name = re.sub("_+", "-", "_".join(s.split("_")[:-1]))
    if "undetermined" in name.lower():
        name = "UND"
    suffix = s.split("_")[-1][0].lower()
    return "_".join([_OWNER_PREFIX, name, suffix])


combined_statistics_df.reset_index(inplace=True)
combined_statistics_df["suggested_strain_name"] = combined_statistics_df[INDEX_COL_NAME].apply(define_strain_name)
combined_statistics_file = os.path.join(ProjectDescriber.SAMPLE_DATA_DIR, "combined_assembly_statistics.tsv")

Utilities.dump_tsv(combined_statistics_df, combined_statistics_file, col_names=[INDEX_COL_NAME] + sorted(
    [i for i in combined_statistics_df.columns if "file" not in i and i not in ("raw_reads", INDEX_COL_NAME)]))

print(combined_statistics_file)
# /data1/bio/projects/vradchenko/lactobacillus_salivarius/sample_data/combined_assembly_statistics.tsv
# Copy the data

ncbi_genome_metadata_df = combined_statistics_df.loc[:, [INDEX_COL_NAME, "real_assembly_coverage",
                                                         "assembly_file"]].copy()
ncbi_genome_metadata_df.rename(columns={"real_assembly_coverage": "genome_coverage"}, inplace=True)

ncbi_genome_metadata_df["assembly_date"] = ncbi_genome_metadata_df["assembly_file"].apply(
    lambda x: datetime.fromtimestamp(os.path.getmtime(x)).strftime("%Y-%m-%d"))
ncbi_genome_metadata_df["assembly_method"] = "SPAdes"
ncbi_genome_metadata_df["sequencing_technology"] = "Illumina MiSeq"
ncbi_genome_metadata_df["filename"] = ncbi_genome_metadata_df["assembly_file"].apply(lambda x: os.path.basename(x))


def parse_spades_version(sample_name_):
    log_file = [i for i in Utilities.scan_whole_dir(
        "/data1/bio/projects/vradchenko/lactobacillus_salivarius/pga-pe/log")
                if i.endswith(".log") and all(j in i for j in ["spades", sample_name_])][0]
    log_lines = Utilities.load_list(log_file)
    image_version_line = [i for i in log_lines if i.strip().startswith("Status: Image is up to date for ")][0].strip()
    spades_version = re.split("[\t ]+", image_version_line)[-1]
    return spades_version


ncbi_genome_metadata_df["assembly_method_version"] = ncbi_genome_metadata_df[INDEX_COL_NAME].apply(
    parse_spades_version)

upload_dir = os.path.join(ProjectDescriber.ROOT_DIR, "assemblies2ncbi")
os.makedirs(upload_dir, exist_ok=True)
ncbi_genome_metadata_df["uploading_assembly_file"] = ncbi_genome_metadata_df["filename"].apply(
    lambda x: os.path.join(upload_dir, x))

ncbi_genome_metadata_file = os.path.join(ProjectDescriber.SAMPLE_DATA_DIR, "sample_batch_genome_accs.tsv")
Utilities.dump_tsv(ncbi_genome_metadata_df, ncbi_genome_metadata_file)
print(ncbi_genome_metadata_file)
# /data1/bio/projects/vradchenko/lactobacillus_salivarius/sample_data/sample_batch_genome_accs.tsv
# Select and upload only required data

source_target_dicts = list(
    ncbi_genome_metadata_df.loc[:, ["assembly_file", "uploading_assembly_file"]].transpose().to_dict().values())
for source_target_dict in source_target_dicts:
    copy2(source_target_dict.get("assembly_file"), source_target_dict.get("uploading_assembly_file"))

print(upload_dir)
# /data1/bio/projects/vradchenko/lactobacillus_salivarius/assemblies2ncbi
# Upload data to NCBI

contamination_reports_dir = os.path.join(upload_dir, "contaminations")
os.makedirs(contamination_reports_dir, exist_ok=True)
# Copy NCBI contamination reports

# Decontaminate assemblies
decontaminated_assemblies_dir = os.path.join(upload_dir, "decontaminated")
os.makedirs(decontaminated_assemblies_dir, exist_ok=True)

for contaminated_assembly in [i.get("uploading_assembly_file") for i in source_target_dicts]:
    contaminated_basename = os.path.splitext(os.path.basename(contaminated_assembly))[0]
    contamination_report = os.path.join(contamination_reports_dir, "Contamination_{}.txt".format(contaminated_basename))
    decontaminated_assembly = os.path.join(decontaminated_assemblies_dir, os.path.basename(contaminated_assembly))
    if not os.path.isfile(contamination_report):
        copy2(contaminated_assembly, decontaminated_assembly)
        continue
    remover = ContaminationRemover(contaminated_assembly, contamination_report)
    remover.export(decontaminated_assembly)

contamination_reports_dir2 = os.path.join(upload_dir, "contaminations2")
os.makedirs(contamination_reports_dir, exist_ok=True)

decontaminated_assemblies_dir2 = os.path.join(upload_dir, "decontaminated2")
os.makedirs(decontaminated_assemblies_dir2, exist_ok=True)

for contaminated_assembly in [i for i in Utilities.scan_whole_dir(decontaminated_assemblies_dir) if i.endswith(".fna")]:
    contaminated_basename = os.path.splitext(os.path.basename(contaminated_assembly))[0]
    contamination_report = os.path.join(contamination_reports_dir2, "Contamination_{}.txt".format(contaminated_basename))
    decontaminated_assembly = os.path.join(decontaminated_assemblies_dir2, os.path.basename(contaminated_assembly))
    if not os.path.isfile(contamination_report):
        continue
    remover = ContaminationRemover(contaminated_assembly, contamination_report)
    remover.export(decontaminated_assembly)
