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
import subprocess
import pandas as pd
from meta.scripts.Utilities import Utilities
from inicolaeva.klebsiella_infants.ProjectDescriber import ProjectDescriber
from Bio import SeqIO
from copy import deepcopy
from meta.scripts.ncbi_contamination_remover import ContaminationRemover
from shutil import copy2

ASSEMBLY_TYPES = ("genome", "plasmid")
ORGANISM = "Klebsiella pneumoniae"
ISOLATE_PREFIX = "KZN_INI_KINF"

assembler_result_dir = "/data1/bio/projects/inicolaeva/klebsiella_infants/pipeline/05_spades"
assembly_files = [i for i in Utilities.scan_whole_dir(assembler_result_dir) if os.path.basename(i) == "contigs.fasta"]
assemblies_target_dir = "/data1/bio/projects/inicolaeva/klebsiella_infants/assemblies"

sample_dirs = sorted(set([os.path.dirname(os.path.dirname(i)) for i in assembly_files]))

_ = subprocess.getoutput("rm -rf {}".format(assemblies_target_dir))
os.makedirs(assemblies_target_dir, exist_ok=True)


assemblies_annotations = []
for sample_dir in sample_dirs:
    sample_name = os.path.basename(sample_dir)
    sample_number = Utilities.safe_findall("([0-9]+)", sample_name)
    sample_assemblies = [i for i in assembly_files if i.startswith(sample_dir)]
    assemblies_annotation = dict()
    seq_records_processed = []
    plasmid_counter = 0
    assembly_target_file = os.path.join(assemblies_target_dir, "{}_genome.fna".format(sample_name))
    for assembly_file_raw in sample_assemblies:
        for assembly_type in ASSEMBLY_TYPES:
            if os.path.dirname(assembly_file_raw).endswith(assembly_type):
                seq_records = sorted(list(SeqIO.parse(assembly_file_raw, "fasta")), key=lambda x: len(x), reverse=True)
                assemblies_annotation["sample_name"] = sample_name
                assemblies_annotation["{}_file".format(assembly_type)] = assembly_file_raw
                seq_records_raw = list(SeqIO.parse(assembly_file_raw, "fasta"))
                assemblies_annotation["{}_assembly_contigs_number_raw".format(assembly_type)] = len(seq_records_raw)
                assemblies_annotation["{}_assembly_bp_raw".format(assembly_type)] = sum(
                    [len(i) for i in seq_records_raw])
                # NCBI does not allow to submit sequences shorter than 200 nucleotides
                seq_records_valid = [i for i in seq_records if len(i) > 200]
                assemblies_annotation["{}_assembly_contigs_number_valid".format(assembly_type)] = len(seq_records_valid)
                assemblies_annotation["{}_assembly_bp_valid".format(assembly_type)] = sum(
                    [len(i) for i in seq_records_valid])
                for seq_record_raw in seq_records_valid:
                    # Example `contigs.fasta` header:
                    # '>NODE_1_length_42950_cov_12.6852_component_0'
                    # contig_number = int(Utilities.safe_findall("^NODE_([0-9]+)", seq_record_raw.id))
                    # Processed FASTA header example:
                    # >contig02 [organism=Clostridium difficile] [strain=ABDC] [plasmid-name=pABDC1] [topology=circular] [completeness=complete]
                    seq_record_processed = deepcopy(seq_record_raw)
                    if assembly_type == "plasmid":
                        plasmid_counter += 1
                        seq_record_processed.description += " PLASMID"
                    seq_records_processed.append(seq_record_processed)
    seq_records_processed = Utilities.remove_duplicate_sequences(seq_records_processed)
    for idx, seq_record_processed in enumerate(seq_records_processed):
        seq_record_processed.id = "contig{a:03d} [organism={b}] [strain={c}_{d}]".format(
            a=idx + 1, b=ORGANISM, c=ISOLATE_PREFIX, d=sample_number)
        if seq_record_processed.description.endswith(" PLASMID"):
            plasmid_counter += 1
            seq_record_processed.description = "[plasmid-name=unnamed{0:02d}]".format(plasmid_counter)
        else:
            seq_record_processed.description = ""
    assemblies_annotations.append(assemblies_annotation)
    #
    SeqIO.write(seq_records_processed, assembly_target_file, "fasta")


INDEX_COL_NAME = "sample_name"
assemblies_statistics_df = pd.DataFrame(assemblies_annotations).set_index(INDEX_COL_NAME)
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

# Decontamination (after NCBI submission)
contamination_reports_dir = os.path.join(assemblies_target_dir, "contamination")
decontaminated_assemblies_dir = os.path.join(assemblies_target_dir, "decontaminated")
_ = [os.makedirs(i, exist_ok=True) for i in (contamination_reports_dir, decontaminated_assemblies_dir)]
# Place reports into this directory

for report_file in Utilities.scan_whole_dir(contamination_reports_dir):
    sample_name_report = Utilities.safe_findall("Contamination_(.+)_genome.txt", report_file)
    for assembly_file in [i for i in Utilities.scan_whole_dir(assemblies_target_dir)
                          if os.path.dirname(i) == assemblies_target_dir and i.endswith("_genome.fna")]:
        sample_name_assembly = Utilities.safe_findall("([^/]+)_genome.fna", assembly_file)
        decontaminated_assembly = os.path.join(
            decontaminated_assemblies_dir, "{}_genome.fna".format(sample_name_assembly))
        if sample_name_report == sample_name_assembly:
            remover = ContaminationRemover(contamination_report=report_file, fna_file=assembly_file)
            remover.export(decontaminated_assembly)
        elif not os.path.isfile(decontaminated_assembly):
            copy2(assembly_file, decontaminated_assembly)

# Add SRA data
ncbi_accessions_df = Utilities.load_tsv("https://raw.githubusercontent.com/ivasilyev/curated_projects/master/inicolaeva/klebsiella_infants/datasets/ncbi_accessions.tsv")

sra_metadata_df = ncbi_accessions_df.loc[:, ["BioSample", "Organism"]]
sra_metadata_df["sample_name"] = "Kleb" + sra_metadata_df["Organism"].str.extract(
    "([0-9]+$)", expand=False)
sra_metadata_df.drop("Organism", axis=1, inplace=True)
sra_metadata_df["title"] = "Illumina-PE-WGS-DNA-Seq of Klebsiella pneumoniae: infant's stool"
sra_metadata_df["library_strategy"] = "WGS"
sra_metadata_df["library_source"] = "GENOMIC"
sra_metadata_df["library_selection"] = "RANDOM"
sra_metadata_df["library_layout"] = "paired"
sra_metadata_df["platform"] = "ILLUMINA"
sra_metadata_df["instrument_model"] = "Illumina MiSeq"
sra_metadata_df["design_description"] = "Libraries were prepared from single colony using the NEBNext Ultra II DNA Library Preparation Kit"
sra_metadata_df["filetype"] = "fastq"

raw_reads_sampledata_df = Utilities.load_tsv(os.path.join(
    os.path.dirname(ProjectDescriber.SAMPLE_DATA_FILE), "raw_reads.sampledata")).set_index(
    "sample_name")
raw_reads_sampledata_df = raw_reads_sampledata_df.loc[:, ["R1", "R2"]].applymap(
    lambda x: os.path.basename(x))
sra_metadata_merged_df = pd.concat([sra_metadata_df.set_index("sample_name"),
                                    raw_reads_sampledata_df], axis=1, sort=False).rename_axis(
    index="library_ID", columns="metadata").reset_index().rename(columns={
        "BioSample": "biosample_accession", "R1": "filename",
        "R2": "filename2"})
sra_metadata_merged_file = os.path.join(os.path.dirname(ProjectDescriber.SAMPLE_DATA_FILE),
                                        "SRA_metadata_kpne.tsv")

# Use the columns order from the corresponding NCBI template
Utilities.dump_tsv(sra_metadata_merged_df, sra_metadata_merged_file, col_names=[i for i in (
    "biosample_accession", "library_ID", "title", "library_strategy", "library_source",
    "library_selection", "library_layout", "platform", "instrument_model", "design_description",
    "filetype", "filename", "filename2", "filename3", "filename4", "assembly", "fasta_file"
) if i in sra_metadata_merged_df.columns])
print(sra_metadata_merged_file)
# /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/SRA_metadata_kpne.tsv
# Copied into the './datasets/ncbi/' directory
