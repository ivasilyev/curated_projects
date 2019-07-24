#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} \
bash

git pull
python3
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from meta.scripts.Utilities import Utilities
from meta.scripts.DigestAssociationsKeeper import DigestAssociationsKeeper
from inicolaeva.klebsiella_infants.ProjectDescriber import ProjectDescriber
from meta.scripts.vfdb.ReferenceDescriber import ReferenceDescriber


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


# Declare class instances & some text constants at global order
projectDescriber = ProjectDescriber()
referenceDescriber = ReferenceDescriber()
digestAssociationsKeeper = DigestAssociationsKeeper()
COVERAGE_DIR = os.path.join(projectDescriber.MAPPED_DATA_DIR, referenceDescriber.ALIAS, "Statistics")
REFERENCE_COL_NAME = "reference_id"
KEYWORDS_ASSOCIATIVE_PAIRS = {"gene_host": {}, "gene_name": digestAssociationsKeeper.VIRULENCE_FACTORS}
DIGEST_LABEL_COL_NAME = "keyword"
RAW_LABEL_COL_NAME = "gene_symbol"
VALUE_COL_NAMES = {"id_mapped_reads_per_million_sample_mapped_reads": "RPM",
                   "id_mapped_reads_per_kbp_per_million_sample_mapped_reads": "RPKM"}
INNER_DONUT_GROUPS = 10
OUTER_DONUT_SUBGROUPS = 5

# Map trimmed reads
# Get VFDB map guideline
get_standalone_map_cmd(projectDescriber.SAMPLE_DATA_FILE, os.path.dirname(COVERAGE_DIR), referenceDescriber)

"""
# Log in into worker node and enter:

export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/nBee.py \
-i /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/trimmed.sampledata \
-r /data/reference/VFDB/vfdb_v2019.04.26/index/vfdb_v2019.04.26_refdata.json \
-o /data1/bio/projects/inicolaeva/klebsiella_infants/mapped/vfdb_v2019.04.26
"""


class InterpretationHandler:
    def __init__(self, reference_describer_instance, value_col_name):
        self.describer = reference_describer_instance
        self.value_col_name = value_col_name
        self.coverage_files = [i for i in Utilities.scan_whole_dir(projectDescriber.MAPPED_DATA_DIR) if
                               all(j in i for j in [self.describer.ALIAS, "coverage.tsv"])]
        self.annotation_file = self.describer.get_refdata_dict().get("sequence_1").annotation_file
        self.raw_annotated_pivot = self.join_and_annotate()
        self.sample_names = self.coverage_files
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
    def update_sample_names(self, regex: str):
        import re
        new_sample_names = [re.findall(regex, i)[0] for i in self.sample_names]
        self.raw_annotated_pivot.rename(columns={i: j for i, j in zip(self.sample_names, new_sample_names)},
                                        inplace=True)
        self.sample_names = new_sample_names
    @staticmethod
    def create_mirrored_path(items: list, sep: str = "_", makedirs: bool = False):
        """
        Simple function to create partially mirroring file path, e.g ["/foo/", "bar", "baz"] -> "/foo/bar/baz/bar_baz"
        :param items: List with root dir as the heading item and other dirs as tailing items
        :param sep: Separator for joining items
        :param makedirs: Should the resulting directory be created?
        :return: file: String with resulting file path (mask)
        """
        file = os.path.join(*items, sep.join(items[1:]))
        if makedirs:
            os.makedirs(os.path.dirname(file), exist_ok=True)
        return file
    @staticmethod
    def make_autopct(values):
        def autopct(pct):
            return "{v}\n({p:.1f}%)".format(v="{:.2g}".format(int(round(pct * sum(values) / 100.0)), "E"), p=pct)
        return autopct


# Create raw coverage pivot for CARD data
for value_col_name in VALUE_COL_NAMES:
    handler = InterpretationHandler(referenceDescriber, value_col_name)
    handler.update_sample_names("\/([A-Z][^\/_]+)_")
    for col_name_with_keywords in KEYWORDS_ASSOCIATIVE_PAIRS:
        df_to_digest = handler.raw_annotated_pivot.loc[:, [col_name_with_keywords] + handler.sample_names]
        associations = KEYWORDS_ASSOCIATIVE_PAIRS.get(col_name_with_keywords)
        if col_name_with_keywords == "gene_host":
            associations = digestAssociationsKeeper.generate_genera_dict(df_to_digest[col_name_with_keywords].values.tolist())
        digest_df, raw_ds = digestAssociationsKeeper.digest_df(df_to_digest, associations=associations, columns_with_keywords=[col_name_with_keywords])
        raw_ds = Utilities.left_merge(raw_ds, handler.raw_annotated_pivot[RAW_LABEL_COL_NAME].reset_index(), REFERENCE_COL_NAME)
        raw_ds[RAW_LABEL_COL_NAME] = raw_ds[RAW_LABEL_COL_NAME].apply(lambda x: min((j for j in str(x).strip().split(" ") if j), key=len))
        for sample_name in digest_df.columns:
            major_digest_df = Utilities.get_n_majors_from_df(digest_df, sample_name, n=INNER_DONUT_GROUPS - 1)
            sample_export_mask = InterpretationHandler.create_mirrored_path([projectDescriber.DATA_DIGEST_DIR, value_col_name, col_name_with_keywords, sample_name], makedirs=True)
            # Create visualization
            fig, ax = plt.subplots()
            _BASE_FONT_SIZE = 15
            plt.rcParams.update({"font.size": _BASE_FONT_SIZE, "figure.figsize": (20, 20)})
            ax.axis("equal")
            y_col_name = major_digest_df.columns[0]
            _WEDGE_WIDTH = 0.3
            _WEDGE_PROPERTIES = dict(width=_WEDGE_WIDTH, edgecolor="w")
            _LABEL_PROPERTIES = dict(fontsize=_BASE_FONT_SIZE, rotation_mode="anchor", verticalalignment="center", horizontalalignment="center")
            # Returning value: [[wedges...], [labels...], [values...]]
            pie_int = ax.pie(major_digest_df[sample_name], radius=1 - _WEDGE_WIDTH, labels=major_digest_df.index,
                             labeldistance=1 - _WEDGE_WIDTH, rotatelabels=False, autopct=InterpretationHandler.make_autopct(major_digest_df[y_col_name]), pctdistance=1 - _WEDGE_WIDTH / 2.0,
                             wedgeprops=_WEDGE_PROPERTIES, textprops=_LABEL_PROPERTIES)
            # Combine color values in 'RGBA' format into the one dictionary
            pie_int_colors = {pie_int[1][idx].get_text(): wedge.get_facecolor() for idx, wedge in enumerate(pie_int[0])}
            # Manual sort the dataset with raw values prior to the order of digest keywords
            major_raw_ds = pd.DataFrame()
            for digest_keyword in major_digest_df.index:
                if digest_keyword == "Other":
                    major_raw_ds_append = pd.DataFrame(major_digest_df.loc["Other"]).transpose()
                    major_raw_ds_append.index.name = DIGEST_LABEL_COL_NAME
                    major_raw_ds_append = major_raw_ds_append.reset_index()
                else:
                    major_raw_ds_append_right = raw_ds.loc[raw_ds[DIGEST_LABEL_COL_NAME] == digest_keyword, [REFERENCE_COL_NAME, sample_name, DIGEST_LABEL_COL_NAME, RAW_LABEL_COL_NAME]]
                    major_raw_ds_append_left = Utilities.get_n_majors_from_df(major_raw_ds_append_right.set_index(REFERENCE_COL_NAME), sample_name, n=OUTER_DONUT_SUBGROUPS - 1).rename(index={"Other": digest_keyword}).reset_index()
                    major_raw_ds_append = Utilities.left_merge(major_raw_ds_append_left, major_raw_ds_append_right, REFERENCE_COL_NAME)
                    major_raw_ds_append[RAW_LABEL_COL_NAME] = major_raw_ds_append[RAW_LABEL_COL_NAME].fillna("{}_Other".format(digest_keyword))
                    major_raw_ds_append[DIGEST_LABEL_COL_NAME] = major_raw_ds_append[DIGEST_LABEL_COL_NAME].fillna("Other")
                pie_ext_append_colors = []
                for row_number in major_raw_ds_append.index.values:
                    row_color = pie_int_colors.get(digest_keyword)
                    if not row_color:
                        continue
                    row_old_alpha = row_color[3]
                    _MINIMAL_ALPHA = 0.2
                    if major_raw_ds_append.shape[0] < 4:
                        row_new_alpha = row_old_alpha - (row_old_alpha * row_number * _MINIMAL_ALPHA)
                    else:
                        row_new_alpha = row_old_alpha - ((row_old_alpha - _MINIMAL_ALPHA) * row_number / float(major_raw_ds_append.shape[0] - 1))
                    pie_ext_append_colors.append(";".join(str(i) for i in list(row_color[:3]) + [row_new_alpha]))
                major_raw_ds_append["color"] = pie_ext_append_colors
                if major_raw_ds_append.shape[0] > 0:
                    if major_raw_ds.shape[0] == 0:
                        major_raw_ds = major_raw_ds_append
                    else:
                        major_raw_ds = pd.concat([major_raw_ds, major_raw_ds_append], axis=0, ignore_index=True, sort=False)
            major_raw_ds = major_raw_ds.fillna("Other")
            pie_ext = ax.pie(major_raw_ds[sample_name], radius=1, labels=major_raw_ds[RAW_LABEL_COL_NAME],
                             labeldistance=1 - _WEDGE_WIDTH / 2, rotatelabels=True, wedgeprops=_WEDGE_PROPERTIES, textprops=_LABEL_PROPERTIES,
                             colors=major_raw_ds["color"].apply(lambda x: tuple(float(i) for i in x.split(";"))).values.tolist())
            # Export visualization tables
            Utilities.dump_tsv(major_digest_df.reset_index(), "{}_inner_values.tsv".format(sample_export_mask))
            Utilities.dump_tsv(major_raw_ds, "{}_outer_values.tsv".format(sample_export_mask))
            # Set labels
            ax.set_xlabel(y_col_name)
            ax.set_ylabel(value_col_name)
            plt.tight_layout()
            # Export PNG
            pie_file = "{}_double_donut.png".format(sample_export_mask)
            fig.suptitle(pie_file, fontsize=_BASE_FONT_SIZE)
            plt.savefig(pie_file, dpi=300, bbox_inches="tight")
            plt.close("all")
            plt.clf()
