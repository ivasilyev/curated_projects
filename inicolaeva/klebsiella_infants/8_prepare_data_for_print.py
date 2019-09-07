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


def process_header(df: pd.DataFrame, capitalize: bool = True):
    _df = df.copy()
    if capitalize:
        _df.rename(columns={i: i[0].upper() + i[1:] for i in _df.columns}, inplace=True)
    return _df.rename(columns={i: i.replace("_", " ") for i in _df.columns})


data_dir = "./inicolaeva/klebsiella_infants/datasets"
article_dir = os.path.join(ProjectDescriber.DATA_DIGEST_DIR, "article")
INDEX_COL_NAME = "sample_name"

antibacterial_agents = ['klebsiella_phage', 'pyo_bacteriophage']
initial_sample_data_df = Utilities.load_tsv(os.path.join(data_dir, "initial_sample_data.tsv")).loc[:, [
    INDEX_COL_NAME] + [
    'sample_number', 'delivery', 'patient_id', 'checkpoint_age_days', 'checkpoint_kpneumoniae_lg_cfu_per_g',
    'extended-spectrum_beta-lactamases'] + antibacterial_agents].set_index(
    INDEX_COL_NAME).sort_index()
initial_sample_data_df["delivery"].replace({"vaginal": "V", "caesarean": "C"}, inplace=True)
initial_sample_data_df["extended-spectrum_beta-lactamases"].replace({True: "+", False: "-"}, inplace=True)
for antibacterial_agent in antibacterial_agents:
    initial_sample_data_df[antibacterial_agent].replace(
        {"susceptible": "S", "resistant": "R", "intermediate": "I", "not defined": "N"}, inplace=True)


initial_sample_data_df.rename(columns={
    "sample_number": "Sample Number", "patient_id": "Patient ID", "checkpoint_age_days": "Patient Age, d",
    "checkpoint_kpneumoniae_lg_cfu_per_g": "K.pneumoniae, lgCFU / g", "extended-spectrum_beta-lactamases": "ESBL",
    'klebsiella_phage': "Klebsiella phage", 'pyo_bacteriophage': "Pyo phage"
}, inplace=True)


antibiotics = sorted([
    'amoxicillin-clavulanic acid', 'ampicillin', 'amikacin', 'aztreonam', 'nitrofurantoin', 'ceftriaxone',
    'sulfamethoxazole', 'trimethoprim', 'ciprofloxacin', 'chloramphenicol', 'fosfomycin', 'netilmicin',
    'gentamicin', 'imipenem', 'meropenem'
])
antibiogram_df = Utilities.load_tsv(
    os.path.join(data_dir, "antibiogram.tsv")).loc[:, [INDEX_COL_NAME] + antibiotics].set_index(INDEX_COL_NAME)
for antibiotic in antibiotics:
    antibiogram_df[antibiotic].replace(
        {"susceptible": "S", "resistant": "R", "intermediate": "I", "not defined": "N"}, inplace=True)


ncbi_accessions_df = Utilities.load_tsv(os.path.join(data_dir, "ncbi_accessions.tsv"))
ncbi_accessions_df[INDEX_COL_NAME] = "Kleb" + ncbi_accessions_df["Organism"].str.extract("(\d+)$")
ncbi_accessions_df = ncbi_accessions_df.loc[:, [INDEX_COL_NAME, "Accession"]].set_index(INDEX_COL_NAME)


combined_assembly_statistics_df = Utilities.load_tsv(
    os.path.join(data_dir, "combined_assembly_statistics.tsv")).loc[:, [INDEX_COL_NAME] + [
        "sample_reads_number", "expected_coverage", "genome_assembly_bp_valid", "genome_assembly_contigs_number_valid"
        ]].set_index(INDEX_COL_NAME)
combined_assembly_statistics_df.rename(columns={
    "sample_reads_number": "Reads Number", "expected_coverage": "Coverage",
    "genome_assembly_bp_valid": "Genome Assembly Length", "genome_assembly_contigs_number_valid": "Contigs Number"
}, inplace=True)


kleborate_results_df = Utilities.load_tsv(
    os.path.join(data_dir, "kleborate_results.tsv")).rename(columns={"strain": INDEX_COL_NAME}).loc[:, [
        INDEX_COL_NAME] + [
        "N50", "largest_contig", "ST", "Yersiniabactin", "Colibactin", "Aerobactin", "Salmochelin", "rmpA", "rmpA2",
        "wzi", "K_locus", "O_locus", "AGly", "Flq", "MLS", "Phe", "Sul", "Tet", "Tmt", "Bla", "Bla_ESBL", "Bla_broad",
        "Bla_broad_inhR"
    ]].set_index(INDEX_COL_NAME)
kleborate_results_df.rename(columns={
    "largest_contig": "Largest Contig Length", "ST": "Sequence Type",
    "genome_assembly_bp_valid": "Genome Assembly Length", "genome_assembly_contigs_number_valid": "Contigs Number",
    "AGly": "Aminoglycosides", "Flq": "Fluoroquinolones", "MLS": "Macrolides", "Phe": "Phenicols",
    "Sul": "Sulfonamides", "Tet": "Tetracyclines", "Tmt": "Trimethoprim", "Bla": "CBL", "Bla_ESBL": "ESBL",
    "Bla_broad": "BSBL", "Bla_broad_inhR": "BSBL-inhR"
}, inplace=True)

phenotype_df = pd.concat([initial_sample_data_df, antibiogram_df], axis=1, sort=False).sort_index()
phenotype_df.index.names = [INDEX_COL_NAME]
phenotype_df = process_header(phenotype_df).transpose().reset_index()
#
Utilities.dump_tsv(phenotype_df, os.path.join(article_dir, "phenotype.tsv"))
Utilities.dump_string(phenotype_df.to_latex(index=False, header=True), os.path.join(article_dir, "phenotype.tex"))


genotype_df = pd.concat([ncbi_accessions_df, combined_assembly_statistics_df, kleborate_results_df],
                        axis=1, sort=False).sort_index()  # .sort_values(["Patient ID", "Sample Number"])
genotype_df.index.names = [INDEX_COL_NAME]
# genotype_df.replace({"_": "\\_"}, regex=True)
genotype_df = process_header(genotype_df, capitalize=False).transpose().reset_index()
#
Utilities.dump_tsv(genotype_df, os.path.join(article_dir, "genotype.tsv"))
Utilities.dump_string(genotype_df.to_latex(index=False, header=True), os.path.join(article_dir, "genotype.tex"))

print(article_dir)
# /data1/bio/projects/inicolaeva/klebsiella_infants/digest/article
