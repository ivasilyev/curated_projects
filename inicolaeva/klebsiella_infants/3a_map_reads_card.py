#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} \
python3
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from meta.scripts.Utilities import Utilities
from meta.scripts.DigestAssociationsKeeper import DigestAssociationsKeeper
from inicolaeva.klebsiella_infants.ProjectDescriber import ProjectDescriber
from meta.scripts.card.ReferenceDescriber import ReferenceDescriber


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
digestAssociationsKeeper = DigestAssociationsKeeper()
COVERAGE_DIR = os.path.join(projectDescriber.MAPPED_DATA_DIR, "Statistics")
REFERENCE_COL_NAME = "reference_id"
VALUE_COL_NAMES = {"id_mapped_reads_per_million_sample_mapped_reads": "RPM",
                   "id_mapped_reads_per_kbp_per_million_sample_mapped_reads": "RPKM"}

# Map trimmed reads
# Get CARD map guideline
get_standalone_map_cmd(projectDescriber.SAMPLE_DATA_FILE, projectDescriber.MAPPED_DATA_DIR, ReferenceDescriber())

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
    def get_gigest_df(self, associations: dict, columns_with_keywords: list):
        association_digest_df = digestAssociationsKeeper.digest_df(
            self.raw_annotated_pivot.loc[:, columns_with_keywords + self.sample_names], associations,
            *columns_with_keywords)
        association_index_name = "_".join(columns_with_keywords)
        association_digest_df.index.name = association_index_name
        association_digest_df.columns.name = "sample_name"
        association_digest_dir = os.path.join(projectDescriber.DATA_DIGEST_DIR, self.describer.ALIAS, value_col_name, association_index_name)
        association_digest_file = os.path.join(association_digest_dir, "digest_{}_{}_{}.tsv".format(self.describer.ALIAS, value_col_name, association_index_name))
        os.makedirs(association_digest_dir, exist_ok=True)
        Utilities.dump_tsv(association_digest_df.reset_index(), table_file=association_digest_file)
        return association_digest_df


# Create raw coverage pivot for CARD data
describer = ReferenceDescriber()
value_col_name = list(VALUE_COL_NAMES)[1]
handler = InterpretationHandler(describer, value_col_name)
genera_names_dict = digestAssociationsKeeper.generate_genera_dict(handler.raw_annotated_pivot["host"].values.tolist())
handler.update_sample_names("\/([A-Z][^\/_]+)_")


columns_with_keywords = ["host"]
df = handler.raw_annotated_pivot.loc[:, columns_with_keywords + handler.sample_names]
associations = digestAssociationsKeeper.generate_genera_dict(df["host"].values.tolist())
digest_df, intermediate_dict = digestAssociationsKeeper.digest_df(df, associations=associations, columns_with_keywords=columns_with_keywords)

major_digest_df = [Utilities.get_n_majors_from_df(digest_df, i) for i in list(digest_df)][-1]

pie_file = os.path.join(projectDescriber.DATA_DIGEST_DIR, "test_pie.png")
os.makedirs(os.path.dirname(pie_file), exist_ok=True)


def make_autopct(values):
    def my_autopct(pct):
        return "{v}\n({p:.1f}%)".format(v="{:.2g}".format(int(round(pct * sum(values) / 100.0)), "E"), p=pct)
    return my_autopct


fig, ax = plt.subplots()
plt.rcParams.update({"font.size": 10, "figure.figsize": (20, 20)})
fig.suptitle(pie_file, fontsize=20)
ax.axis("equal")
y_col_name = major_digest_df.columns[0]

raw_values_dfs_dict = {}
raw_values_dict = {"x": [], "y": []}
for keyword in major_digest_df.index:
    _keyword_pivot_df = intermediate_dict.get(keyword)
    if _keyword_pivot_df is not None:
        # keyword_pivot_df = keyword_pivot_df.loc[keyword_pivot_df[y_col_name].astype(float) > 0.0][y_col_name]
        keyword_pivot_df = Utilities.get_n_majors_from_df(_keyword_pivot_df, col_name=y_col_name, n=5)
        if keyword_pivot_df.shape[0] > 0:
            keyword_pivot_df = keyword_pivot_df.join(handler.raw_annotated_pivot["gene"], how="left").fillna("Other")
            raw_values_dfs_dict[keyword] = keyword_pivot_df
            raw_values_dict["x"].extend([min((j for j in i.strip().split(" ") if j), key=len) for i in keyword_pivot_df["gene"].values.tolist()])
            raw_values_dict["y"].extend(keyword_pivot_df[y_col_name])


size = 0.3
wedgeprops = dict(width=size, edgecolor="w")
pie_in = ax.pie(major_digest_df[y_col_name], radius=1 - size, labels=major_digest_df.index, labeldistance=1 - size,
                autopct=make_autopct(major_digest_df[y_col_name]), wedgeprops=wedgeprops)
pie_out = ax.pie(raw_values_dict.get("y"), radius=1, labels=raw_values_dict.get("x"),
                 wedgeprops=wedgeprops)
ax.set_xlabel(y_col_name)
ax.set_ylabel(value_col_name)
plt.tight_layout()
plt.savefig(pie_file, dpi=300, bbox_inches="tight")
plt.close("all")
plt.clf()
















def make_autopct(values):
    def my_autopct(pct):
        return "{v}\n({p:.1f}%)".format(v="{:.2g}".format(int(round(pct * sum(values) / 100.0)), "E"), p=pct)
    return my_autopct


fig = plt.figure()
sns.set(style="whitegrid", font_scale=1)
y_col_name = major_digest_df.columns[0]

ax = major_digest_df.plot.pie(y=y_col_name, figsize=(10, 10), title=pie_file, autopct=make_autopct(major_digest_df[y_col_name]))
ax.set_xlabel(y_col_name)
ax.set_ylabel(value_col_name)
plt.axis("equal")
plt.tight_layout()
plt.savefig(pie_file, dpi=300, bbox_inches="tight")
plt.close("all")
plt.clf()



