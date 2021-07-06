#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# To access a JupyterLab server on http://ip_address:61156/?token=TOKEN:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -p 61156:61156 -it ${IMG} bash

git pull && jupyter lab --ip=0.0.0.0 --port=61156 --no-browser --NotebookApp.token=TOKEN
"""

import os
import xlrd
import joblib
import numpy as np
import pandas as pd
from itertools import product, combinations_with_replacement
from meta.scripts.Utilities import Utilities
from meta.scripts.utils.diversity_utils import count_alpha_diversity
from meta.scripts.utils.pandas_utils import load_tsv, dump_tsv, concat
from ashestopalov.nutrition.obesity_metagenomes.ProjectDescriber import ProjectDescriber


def mp_apply_function_to_df(func, df: pd.DataFrame):
    results = joblib.Parallel(n_jobs=-1)(joblib.delayed(func)(j) for j in [df[i] for i in df.columns])
    return pd.DataFrame(results)


def add_prefix_to_columns(prefix: str, df: pd.DataFrame):
    return df.rename(columns={i: "{}@{}".format(prefix, i) for i in df.columns})


def select_data_columns(df: pd.DataFrame):
    return df.loc[:, [i for i in df.columns if "@" in i]]


SN = "sample_name"
TA = "target_alias"
source_raw_reads_data_df = load_tsv(os.path.join(ProjectDescriber.SAMPLE_DATA_DIR, "source_raw_reads_data.tsv"))

thinned_sample_names = []
for sample_name in source_raw_reads_data_df[SN].unique():
    sample_name_sub_df = source_raw_reads_data_df.loc[source_raw_reads_data_df[SN] == sample_name]
    top_alias_name = sample_name_sub_df.loc[:, ["target_alias", "file_size"]].groupby(["target_alias"]).sum().sort_values("file_size", ascending=False).head(1).index[0]
    thinned_sample_names.append(dict(sample_name=sample_name, target_alias=top_alias_name))

cured_group_data_df = load_tsv(os.path.join(ProjectDescriber.SAMPLE_DATA_DIR, "cured_group_data.tsv"))
combined_group_data_df = pd.concat([i.set_index(SN) for i in (cured_group_data_df, pd.DataFrame(thinned_sample_names))], sort=False, axis=1).rename_axis(index=SN).dropna()

dump_tsv(combined_group_data_df, os.path.join(ProjectDescriber.SAMPLE_DATA_DIR, "combined_group_data.tsv"), reset_index=True)

AGE_GROUPS = ("adult", "child")
DIAGNOSIS_GROUPS = ("normal", "obesity")
sample_source = "stool"

ASSOCIATIONS = {
    "Indoles": "Indol	Quinolinic	Kynurenine	3HIAA	Antranillic	Xanturenic	Kynurenic	Indol-3-Lactic	Indole-3-acetic	Indol-carboxaldehyde	Indol-acrilyc	Indol-propionic	Tryptamine",
    "ELISA": "VEGF	Adiponectin	Resistin	ASPR	FGF21	Irisin	Myostatin/GDF8	HOMA	Insulin	Leptin",
    "Anthropometry": "Body Mass Index	Waist Circumference",
    "Blood analysis": "Blood Glucose	High-Density Lipoproteins	Thyroglobulin	Total Blood Cholesterol	Low-Density Lipoproteins	Systolic Blood Pressure	Diastolic Blood Pressure",
    "Resorcinols": "Resorcinol	Methylresorcinol	Ethylresorcinol	Propylresorcinol	Penthylresorcinol	Hexylresorcinol	Dodecylresorcinol	Pentadecylresorcinol"
}

associations = {k: v.split("\t") for k, v in ASSOCIATIONS.items()}
associations_to_rename = dict()
for association in associations:
    for word in associations[association]:
        associations_to_rename[word] = "{}@{}".format(association, word)

supplied_data_dir = os.path.join(ProjectDescriber.DATA_DIR, "supplied_data")

raw_data_dfs = dict()
for analysis_type in ("Indoles", "Resorcinols"):
    raw_data_dfs[analysis_type] = pd.read_excel(
        os.path.join(supplied_data_dir, "cured_data_{}.xlsx".format(analysis_type.lower())), sheet_name=sample_source)
    raw_data_dfs[analysis_type].rename(columns=associations_to_rename, inplace=True)

merged_data_dir = os.path.join(ProjectDescriber.DATA_DIR, "merged_data")

feature_dfs = dict()
for feature_name in ("EC", "KO", "OTU", "pathway"):
    feature_df = load_tsv(os.path.join(merged_data_dir, "{}_IDs.tsv".format(feature_name)))
    id_column_names = np.intersect1d(feature_df.columns,
                                     ("description", "function", "#OTU ID", "taxonomy", "pathway"))
    feature_df_index_name = "_".join(id_column_names)
    feature_df[feature_df_index_name] = feature_df.loc[:, id_column_names].apply(
        lambda x: "&".join(["{}={}".format(x.index[idx], i) for idx, i in enumerate(x.astype(str))]), axis=1)
    feature_df = feature_df.drop(id_column_names, axis=1).set_index(feature_df_index_name)
    feature_dfs[feature_name] = add_prefix_to_columns(
        feature_name, feature_df.transpose().rename_axis(index=TA))

otu_df = load_tsv(os.path.join(merged_data_dir, "OTU_IDs.tsv"))
# Setting taxonomy as index in order to count the Faith estimators
otu_samples_df = otu_df.set_index("taxonomy").drop("#OTU ID", axis=1)
alpha_diversity_df = mp_apply_function_to_df(count_alpha_diversity, otu_samples_df)
alpha_diversity_df = add_prefix_to_columns("Alpha", alpha_diversity_df)
feature_dfs["Alpha"] = alpha_diversity_df

for feature_name, feature_df in feature_dfs.items():
    feature_df = concat([combined_group_data_df.set_index("target_alias", drop=False), feature_df])
    feature_dfs[feature_name] = feature_df.query("sample_source == '{}'".format(sample_source))

feature_dfs.update(raw_data_dfs)

correlation_dir = os.path.join(ProjectDescriber.DATA_DIR, "correlation_data", "group_datasets")

group_combinations = list(product(combinations_with_replacement(feature_dfs.keys(), 2), AGE_GROUPS, DIAGNOSIS_GROUPS))

correlation_tables = []
for feature_pair, age, diagnosis in group_combinations:
    query = "age == '{}' and diagnosis == '{}'".format(age, diagnosis)
    correlation_df = concat([select_data_columns(feature_dfs[i].query(query).reset_index()) for i in sorted(set(feature_pair))])
    correlation_table = os.path.join(correlation_dir, "{}_for_{}_{}.tsv".format(" vs ".join(feature_pair), age, diagnosis))
    dump_tsv(correlation_df, correlation_table)
    correlation_tables.append(correlation_table)

Utilities.dump_list(correlation_tables, os.path.join(correlation_dir, "tables.txt"))
Utilities.dump_list(correlation_tables, os.path.join(correlation_dir, "tables.txt.bak"))
