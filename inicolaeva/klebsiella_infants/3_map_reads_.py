#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from meta.scripts.Utilities import Utilities
from meta.scripts.DigestAssociationsKeeper import DigestAssociationsKeeper
from inicolaeva.klebsiella_infants.ProjectDescriber import ProjectDescriber
from meta.scripts.card.ReferenceDescriber import ReferenceDescriber as cardDescriber


def get_standalone_map_cmd(input_sampledata: str, output_dir: str, describer_instance):
    print("""
# Log in into worker node and enter:

export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \\
docker pull $IMG && \\
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \\
python3 /home/docker/scripts/nBee.py \\
-i {sd} \\
-r {rd} \\
-o {out}
""".format(sd=input_sampledata, rd=describer_instance.REFDATA, out=output_dir))


# Declare class instances & some text constants
projectDescriber = ProjectDescriber()
associationsKeeper = DigestAssociationsKeeper()
COVERAGE_DIR = os.path.join(projectDescriber.MAPPED_DATA_DIR, "Statistics")
REFERENCE_COL_NAME = "reference_id"
VALUE_COL_NAMES = {"id_mapped_reads_per_million_sample_mapped_reads": "RPM",
                   "id_mapped_reads_per_kbp_per_million_sample_mapped_reads": "RPKM"}

# Map trimmed reads
# Get CARD map guideline
get_standalone_map_cmd(projectDescriber.SAMPLE_DATA_FILE, projectDescriber.MAPPED_DATA_DIR, cardDescriber())

"""
# Log in into worker node and enter:

export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/nBee.py \
-i /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/trimmed.sampledata \
-r /data/reference/CARD/card_v3.0.1/index/card_v3.0.1_refdata.json \
-o /data1/bio/projects/inicolaeva/klebsiella_infants/mapped
"""


class InterpretationHandler:
    def __init__(self, reference_describer_instance, value_col_name):
        self.describer = reference_describer_instance
        self.value_col_name = value_col_name
        self.coverage_files = [i for i in Utilities.scan_whole_dir(projectDescriber.MAPPED_DATA_DIR) if
                               all(j in i for j in [self.describer.ALIAS, "coverage.tsv"])]
        self.annotation_file = self.describer.get_refdata_dict().get("sequence_1").annotation_file
        self.raw_annotated_pivot = self.join_and_annotate()
    @staticmethod
    def join_by_value_columns(tables: list, index_col_name: str, value_col_name_: str):
        dfs_list = [Utilities.load_tsv(i).set_index(index_col_name)[value_col_name_].rename(i) for i in tables]
        out = pd.concat(dfs_list, axis=1, sort=False).sort_index()
        out.index.names = [index_col_name]
        return out
    def join_and_annotate(self):
        annotation_df = Utilities.load_tsv(self.annotation_file).set_index(REFERENCE_COL_NAME)
        values_df = self.join_by_value_columns(self.coverage_files, REFERENCE_COL_NAME, self.value_col_name)
        out = pd.concat([annotation_df, values_df], axis=1, sort=False)
        out.index.names = [REFERENCE_COL_NAME]
        return out
    @staticmethod
    def get_genera_dict(input_list: list):
        return {j: () for j in sorted(
            [Utilities.safe_findall("([A-Z][a-z]{4,})", i).strip() for i in set(input_list) if isinstance(i, str)]) if
                len(j) > 0}








# Create raw coverage pivot for CARD data
card_describer = cardDescriber()
value_col_name = list(VALUE_COL_NAMES)[1]
card_handler = InterpretationHandler(card_describer, value_col_name)
genera_names_dict = card_handler.get_genera_dict(card_handler.raw_annotated_pivot["host"].values.tolist())



import os
import re
import matplotlib as mpl
from meta.scripts.Utilities import Utilities
from meta.scripts.DigestAssociationsKeeper import DigestAssociationsKeeper

PROJECT_ROOT_DIR = "/data1/bio/projects/inicolaeva/klebsiella_infants"
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
-i /data1/bio/projects/inicolaeva/klebsiella_infants/trimmed.sampledata \
-r /data/reference/CARD/card_v3.0.1/index/card_v3.0.1_refdata.json \
-o /data1/bio/projects/inicolaeva/klebsiella_infants/map_data
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
        association_digest.reset_index().to_csv(
            os.path.join(card_digest_dir, "digest_card_{}_{}.tsv".format(value_col_name, annotation_col_name)),
            sep="\t", index=False, header=True)
        association_digest_percentage = association_digest * 100 / association_digest.sum()
        fig = plt.figure()
        sns.set(style="whitegrid", font_scale=1)
        export_df = association_digest.rename(columns={i: Utilities.safe_findall(
            "/data1/bio/projects/inicolaeva/klebsiella_infants/map_data/Statistics/(.+)_card_v3.0.1_coverage.tsv", i)
        for i in
            list(association_digest)}).transpose()
        export_df.index.name = "sample_name"
        ax = export_df.plot(kind='bar', stacked='True', figsize=(20, 10))
        ax.set_ylabel(value_col_name)
        legend = ax.legend(loc="center left", shadow=True, fontsize="x-small", bbox_to_anchor=(1.04, 0.5),
                           borderaxespad=0)
        image_file_name = os.path.join(card_digest_dir,
                                       "digest_card_{}_{}.png".format(value_col_name, annotation_col_name))
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
-i /data1/bio/projects/inicolaeva/klebsiella_infants/trimmed.sampledata \
-r /data/reference/TADB/tadb_v2.0/index/tadb_v2.0_refdata.json \
-o /data1/bio/projects/inicolaeva/klebsiella_infants/map_data
"""

tadb_coverage_files = [i for i in Utilities.scan_whole_dir(MAPPED_STAT_DIR) if all(j in i for j in ["tadb", "coverage.tsv"])]

for value_col_name in VALUE_COL_NAMES:
    tadb_annotated_df = join_and_annotate(tadb_coverage_files, value_col_name,
                                          "/data/reference/TADB/tadb_v2.0/index/tadb_v2.0_annotation.tsv")
    tadb_raw_dir = os.path.join(RAW_DIR, "tadb")
    os.makedirs(tadb_raw_dir, exist_ok=True)
    tadb_annotated_df.reset_index().to_csv(
        os.path.join(tadb_raw_dir, "tadb_annotated_pivot_by_{}.tsv".format(value_col_name)), sep="\t", index=False,
        header=True)
    genera_names_dict = generate_keywords_dict(
        tadb_annotated_df["host"].str.extract("([A-Z][a-z]+)")[0].values.tolist())
    for annotation_col_name, keywords_dict in zip(["protein_description", "host"],
                                                  [keeper.VIRULENCE_FACTORS, genera_names_dict]):
        association_digest = keeper.digest_df(tadb_annotated_df.loc[:, tadb_coverage_files + [annotation_col_name]],
                                              keywords_dict, annotation_col_name)
        tadb_digest_dir = os.path.join(DIGEST_DIR, "tadb", value_col_name, annotation_col_name)
        os.makedirs(tadb_digest_dir, exist_ok=True)
        association_digest.reset_index().to_csv(
            os.path.join(tadb_digest_dir, "digest_tadb_{}_{}.tsv".format(value_col_name, annotation_col_name)),
            sep="\t", index=False, header=True)
        association_digest_percentage = association_digest * 100 / association_digest.sum()
        fig = plt.figure()
        sns.set(style="whitegrid", font_scale=1)
        export_df = association_digest.rename(columns={i: Utilities.safe_findall(
            "/data1/bio/projects/inicolaeva/klebsiella_infants/map_data/Statistics/(.+)_tadb_v2.0_coverage.tsv", i) for
        i in
            list(association_digest)}).transpose()
        export_df.index.name = "sample_name"
        ax = export_df.plot(kind='bar', stacked='True', figsize=(20, 10))
        ax.set_ylabel(value_col_name)
        legend = ax.legend(loc="center left", shadow=True, fontsize="x-small", bbox_to_anchor=(1.04, 0.5),
                           borderaxespad=0)
        image_file_name = os.path.join(tadb_digest_dir,
                                       "digest_tadb_{}_{}.png".format(value_col_name, annotation_col_name))
        ax.set_title(os.path.basename(image_file_name), fontsize="large")
        plt.savefig(image_file_name, dpi=300, bbox_inches="tight")
        plt.close()
        plt.clf()

    @staticmethod
    def create_mirrored_path(items: list, sep: str = "_", create_dirs: bool = False):
        """
        Simple function to create partially mirroring file path, e.g ["/foo/", "bar", "baz"] -> "/foo/bar/baz/bar_baz"
        :param items: List with root dir as the heading item and other dirs as tailing items
        :param sep: Separator for joining items
        :param create_dirs: Should the resulting directory be created?
        :return: file: String with resulting file path (mask)
        """
        file = os.path.join(*items, sep.join(items[1:]))
        if create_dirs:
            os.makedirs(os.path.dirname(file), exist_ok=True)
        return file
