#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3
"""

import os
import re
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from meta.scripts.Utilities import Utilities
from meta.scripts.DigestAssociationsKeeper import DigestAssociationsKeeper

PROJECT_ROOT_DIR = "/data1/bio/projects/auhrbach/klebsiella_infants"
MAPPED_STAT_DIR = os.path.join(PROJECT_ROOT_DIR, "map_data", "Statistics")
RAW_DIR = os.path.join(PROJECT_ROOT_DIR, "raw")
DIGEST_DIR = os.path.join(PROJECT_ROOT_DIR, "digest")

REFERENCE_COL_NAME = "reference_id"
VALUE_COL_NAMES = ("id_mapped_reads_per_million_sample_mapped_reads",
                   "id_mapped_reads_per_kbp_per_million_sample_mapped_reads")
keeper = DigestAssociationsKeeper()


# Map CARD data
"""
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/nBee.py \
-i /data1/bio/projects/auhrbach/klebsiella_infants/trimmed.sampledata \
-r /data/reference/CARD/card_v3.0.1/index/card_v3.0.1_refdata.json \
-o /data1/bio/projects/auhrbach/klebsiella_infants/map_data
"""


def read_table(table_file: str, index=None):
    df = pd.read_table(table_file, header=0, sep="\t", encoding="utf-8")
    if index:
        df = df.set_index(index)
    return df


def join_by_value_columns(tables: list, index_col_name: str, _value_col_name: str):
    dfs_list = [read_table(i, index_col_name)[_value_col_name].rename(i) for i in tables]
    return pd.concat(dfs_list, axis=1, sort=False).sort_index()


def join_and_annotate(table_files: list, _value_col_name, annotation_file: str):
    annotation_df = read_table(annotation_file, REFERENCE_COL_NAME)
    values_df = join_by_value_columns(table_files, REFERENCE_COL_NAME, _value_col_name)
    return pd.concat([annotation_df, values_df], axis=1, sort=False)


def generate_keywords_dict(keywords: list):
    return {j: () for j in sorted([i for i in set(keywords) if isinstance(i, str)])}


card_coverage_files = [i for i in Utilities.ls(MAPPED_STAT_DIR) if all(j in i for j in ["card", "coverage.tsv"])]

for value_col_name in VALUE_COL_NAMES:
    card_annotated_df = join_and_annotate(card_coverage_files, value_col_name,
                                          "/data/reference/CARD/card_v3.0.1/index/card_v3.0.1_annotation.tsv")
    card_raw_dir = os.path.join(RAW_DIR, "card")
    os.makedirs(card_raw_dir, exist_ok=True)
    card_annotated_df.reset_index().to_csv(
        os.path.join(card_raw_dir, "card_annotated_pivot_by_{}.tsv".format(value_col_name)), sep="\t", index=False,
        header=True)
    genera_names_dict = {i: () for i in sorted(
        [i for i in set(card_annotated_df["host"].str.extract("([A-Z][a-z]+)")[0].values.tolist()) if
         isinstance(i, str)])}
    for annotation_col_name, keywords_dict in zip(["Resistance Mechanism", "Drug Class", "host"],
                                                  [keeper.RESISTANCE_MECHANISMS, keeper.DRUG_CLASSES,
                                                   genera_names_dict]):
        association_digest = keeper.digest_df(card_annotated_df.loc[:, card_coverage_files + [annotation_col_name]],
                                              keywords_dict, annotation_col_name)
        card_digest_dir = os.path.join(DIGEST_DIR, "card", value_col_name, annotation_col_name)
        os.makedirs(card_digest_dir, exist_ok=True)
        association_digest.reset_index().to_csv(os.path.join(card_digest_dir, "digest_card_{}_{}.tsv".format(value_col_name, annotation_col_name)), sep="\t", index=False, header=True)
        association_digest_percentage = association_digest * 100 / association_digest.sum()
        fig = plt.figure()
        sns.set(style="whitegrid", font_scale=1)
        export_df = association_digest.rename(columns={i: Utilities.safe_findall(
            "/data1/bio/projects/auhrbach/klebsiella_infants/map_data/Statistics/(.+)_card_v3.0.1_coverage.tsv", i) for i in
                                                       list(association_digest)}).transpose()
        export_df.index.name = "sample_name"
        ax = export_df.plot(kind='bar', stacked='True', figsize=(20, 10))
        ax.set_ylabel(value_col_name)
        legend = ax.legend(loc="center left", shadow=True, fontsize="x-small", bbox_to_anchor=(1.04, 0.5), borderaxespad=0)
        image_file_name = os.path.join(card_digest_dir, "digest_card_{}_{}.png".format(value_col_name, annotation_col_name))
        ax.set_title(os.path.basename(image_file_name), fontsize="large")
        plt.savefig(image_file_name, dpi=300, bbox_inches="tight")
        plt.close()
        plt.clf()


# Map TADB data
"""
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/nBee.py \
-i /data1/bio/projects/auhrbach/klebsiella_infants/trimmed.sampledata \
-r /data/reference/TADB/tadb_v2.0/index/tadb_v2.0_refdata.json \
-o /data1/bio/projects/auhrbach/klebsiella_infants/map_data
"""

tadb_coverage_files = [i for i in Utilities.ls(MAPPED_STAT_DIR) if all(j in i for j in ["tadb", "coverage.tsv"])]

for value_col_name in VALUE_COL_NAMES:
    tadb_annotated_df = join_and_annotate(tadb_coverage_files, value_col_name,
                                          "/data/reference/TADB/tadb_v2.0/index/tadb_v2.0_annotation.tsv")
    tadb_raw_dir = os.path.join(RAW_DIR, "tadb")
    os.makedirs(tadb_raw_dir, exist_ok=True)
    tadb_annotated_df.reset_index().to_csv(
        os.path.join(tadb_raw_dir, "tadb_annotated_pivot_by_{}.tsv".format(value_col_name)), sep="\t", index=False,
        header=True)
    genera_names_dict = generate_keywords_dict(tadb_annotated_df["host"].str.extract("([A-Z][a-z]+)")[0].values.tolist())
    for annotation_col_name, keywords_dict in zip(["protein_description", "host"],
                                                  [keeper.VIRULENCE_FACTORS, genera_names_dict]):
        association_digest = keeper.digest_df(tadb_annotated_df.loc[:, tadb_coverage_files + [annotation_col_name]],
                                              keywords_dict, annotation_col_name)
        tadb_digest_dir = os.path.join(DIGEST_DIR, "tadb", value_col_name, annotation_col_name)
        os.makedirs(tadb_digest_dir, exist_ok=True)
        association_digest.reset_index().to_csv(os.path.join(tadb_digest_dir, "digest_tadb_{}_{}.tsv".format(value_col_name, annotation_col_name)), sep="\t", index=False, header=True)
        association_digest_percentage = association_digest * 100 / association_digest.sum()
        fig = plt.figure()
        sns.set(style="whitegrid", font_scale=1)
        export_df = association_digest.rename(columns={i: Utilities.safe_findall(
            "/data1/bio/projects/auhrbach/klebsiella_infants/map_data/Statistics/(.+)_tadb_v2.0_coverage.tsv", i) for i in
                                                       list(association_digest)}).transpose()
        export_df.index.name = "sample_name"
        ax = export_df.plot(kind='bar', stacked='True', figsize=(20, 10))
        ax.set_ylabel(value_col_name)
        legend = ax.legend(loc="center left", shadow=True, fontsize="x-small", bbox_to_anchor=(1.04, 0.5), borderaxespad=0)
        image_file_name = os.path.join(tadb_digest_dir, "digest_tadb_{}_{}.png".format(value_col_name, annotation_col_name))
        ax.set_title(os.path.basename(image_file_name), fontsize="large")
        plt.savefig(image_file_name, dpi=300, bbox_inches="tight")
        plt.close()
        plt.clf()
