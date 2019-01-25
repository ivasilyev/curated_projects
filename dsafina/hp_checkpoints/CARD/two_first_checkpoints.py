#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3
"""

from meta.scripts.GroupDataPreparer import GroupDataAssemblyGuideLiner
from meta.scripts.DigestAssociationsKeeper import DigestAssociationsKeeper
from meta.scripts.PivotSplitter import PivotSplitter
import subprocess
import os
import re
from meta.scripts.Utilities import Utilities
from dsafina.hp_checkpoints.ProjectDescriber import ProjectDescriber
from meta.scripts.card.ReferenceDescriber import ReferenceDescriber
from meta.scripts.LaunchGuideLiner import LaunchGuideLiner
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import xlrd

projectDescriber = ProjectDescriber()

requestedSamplePairs = """
99-100
106-107
108-109
112-113
114-115
116-117
118-119
120-121
122-123
124-125
126-127
128-129
130-131
132-133
134-135
188-189
190-191
50-51
55-56
57-58
59-60
61-62
67-68
71-72
75-76
79-80
89-90
91-92
97-98
"""

requestedContolGroupSamples = [12, 14, 154, 166, 17, 20, 30, 36, 37, 44, 45, 46] + list(range(141, 152)) + list(range(174, 184)) + list(range(203, 212))
requestedContolGroupDF = pd.DataFrame([{"sample": "{}HP".format(i), "group": "control"} for i in requestedContolGroupSamples]).loc[:, ["sample", "group"]]

requestedGroupDataDF = pd.DataFrame([{"1st": l[0], "2nd": l[1]} for l in [["{}HP".format(k) for k in j.split("-")] for j in [i.strip() for i in re.sub("\r\n+", "\n", requestedSamplePairs).split("\n") if len(i.strip()) > 0]]]).melt().rename(columns={"value": "sample", "variable": "group"}).loc[:, ["sample", "group"]]
requestedGroupDataDF = pd.concat([requestedGroupDataDF, requestedContolGroupDF], axis=0, ignore_index=True)
group_names = sorted(set(requestedGroupDataDF["group"]))

projectDescriber.groupdata = "/data1/bio/projects/dsafina/hp_checkpoints/two_first_checkpoints.groupdata"
requestedGroupDataDF.to_csv(projectDescriber.groupdata, sep="\t", index=False, header=False)

index_col_name = "reference_id"
outputDir = "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/two_first_checkpoints/"
prefix = "/data2/bio/Metagenomes/CARD/Statistics/"
suffix = "_card_v2.0.3_coverage.tsv"

# RPM and RPKM
value_col_names = ("id_mapped_reads_per_million_sample_mapped_reads",
                   "id_mapped_reads_per_kbp_per_million_sample_mapped_reads")
digest_value_col_names = ("log2(RPM+1)", "kRPKM")

for value_col_name in value_col_names:
    guideliner = GroupDataAssemblyGuideLiner(groupdata=projectDescriber.groupdata,
                                             index_column=index_col_name,
                                             output_dir="{}single_genes/{}/".format(outputDir, value_col_name),
                                             prefix=prefix,
                                             suffix=suffix,
                                             value_column=value_col_name)
    print("{}\n".format(guideliner.external_launch_command))

subprocess.getoutput("rm -rf {}".format(outputDir))

"""
# Pre-setup to launch from different node for group data file "/data1/bio/projects/dsafina/hp_checkpoints/two_first_checkpoints.groupdata" and value column "id_mapped_reads_per_million_sample_mapped_reads":

export IMG=ivasilyev/curated_projects:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it $IMG python3 /home/docker/scripts/statistical_tools/groupdata2statistics.py \
-g "/data1/bio/projects/dsafina/hp_checkpoints/two_first_checkpoints.groupdata" \
-p "/data2/bio/Metagenomes/CARD/Statistics/" \
-s "_card_v2.0.3_coverage.tsv" \
-i "reference_id" \
-v "id_mapped_reads_per_million_sample_mapped_reads" \
-o "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/two_first_checkpoints/single_genes/id_mapped_reads_per_million_sample_mapped_reads/"



# Pre-setup to launch from different node for group data file "/data1/bio/projects/dsafina/hp_checkpoints/two_first_checkpoints.groupdata" and value column "id_mapped_reads_per_kbp_per_million_sample_mapped_reads":

export IMG=ivasilyev/curated_projects:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it $IMG python3 /home/docker/scripts/statistical_tools/groupdata2statistics.py \
-g "/data1/bio/projects/dsafina/hp_checkpoints/two_first_checkpoints.groupdata" \
-p "/data2/bio/Metagenomes/CARD/Statistics/" \
-s "_card_v2.0.3_coverage.tsv" \
-i "reference_id" \
-v "id_mapped_reads_per_kbp_per_million_sample_mapped_reads" \
-o "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/two_first_checkpoints/single_genes/id_mapped_reads_per_kbp_per_million_sample_mapped_reads/"
"""

annotation_df = pd.read_table("/data/reference/CARD/card_v2.0.3/index/card_v2.0.3_annotation.tsv").set_index(index_col_name)

subprocess.getoutput("rm -rf /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/two_first_checkpoints/gene_groups")


class DigestKeeper:
    def __init__(self):
        self.value_col_name = ""
        self.digest_value_col_name = ""
        self.samples_dir = ""
        self.groupdata_df = pd.DataFrame()
        self.single_genes_ds = pd.DataFrame()
        self.gene_groups_ds = pd.DataFrame()


digest_keepers_dict = {}
for value_col_name, digest_value_col_name in zip(value_col_names, digest_value_col_names):
    digest_keeper = DigestKeeper()
    digest_keeper.value_col_name, digest_keeper.digest_value_col_name = value_col_name, digest_value_col_name
    digest_keeper.samples_dir = "{a}gene_groups/{b}/".format(a=outputDir, b=digest_value_col_name)
    os.makedirs(digest_keeper.samples_dir, exist_ok=True)
    for group_name in group_names:
        single_genes_pivot_file = "{a}single_genes/{b}/pivot_by_{b}_{c}.tsv".format(a=outputDir, b=value_col_name, c=group_name)
        single_genes_pivot_df = pd.read_table(single_genes_pivot_file).set_index(index_col_name)
        digest_keeper.single_genes_df = pd.concat([annotation_df, single_genes_pivot_df], axis=1)
        # Annotate pivot
        single_genes_annotated_df = pd.concat([annotation_df, single_genes_pivot_df], axis=1)
        # Name the df
        single_genes_annotated_df.name, single_genes_annotated_df.columns.name = value_col_name, "sample_path"
        # Sum all rows with keywords (dict values) in proper column and assign a single name (dict key) to the sum row
        gene_groups_total_df = pd.DataFrame()
        for keywords_col_name, association in zip(("Drug Class", "Resistance Mechanism"), (DigestAssociationsKeeper.DRUG_CLASSES, DigestAssociationsKeeper.RESISTANCE_MECHANISMS)):
            single_genes_keywords_df = digest_keeper.single_genes_df.loc[:, [keywords_col_name] + list(single_genes_pivot_df)]
            gene_groups_keywords_df = DigestAssociationsKeeper.digest_df(df=single_genes_keywords_df, associations=association, columns_with_keywords=[keywords_col_name])
            # Prepare values for visualization
            if digest_value_col_name == "log2(RPM+1)":
                gene_groups_keywords_df = gene_groups_keywords_df.apply(lambda x: np.log2(x + 1))
            elif digest_value_col_name == "kRPKM":
                gene_groups_keywords_df = gene_groups_keywords_df.apply(lambda x: x.astype(float) / 1000.0)
            # Sort each keyword group and join vertically
            gene_groups_total_df = pd.concat([gene_groups_total_df, gene_groups_keywords_df.sort_index()], axis=0)
        # Name the df
        gene_groups_total_df.name, gene_groups_total_df.columns.name = digest_value_col_name, "sample_path"
        splitter = PivotSplitter(pivot_df=gene_groups_total_df, value_col_name=digest_value_col_name)
        splitter.split(output_dir="{}samples/{}/".format(digest_keeper.samples_dir, group_name))
        digest_keeper.groupdata_df = pd.concat([digest_keeper.groupdata_df, splitter.get_groupdata(group_name)], axis=0, ignore_index=True)
        # Prepare datasets
        single_genes_total_ds = single_genes_pivot_df.reset_index().melt(id_vars=[index_col_name], var_name="sample_path").rename(columns={"value": value_col_name})
        single_genes_total_ds["group_name"] = group_name
        single_genes_total_ds = single_genes_total_ds.set_index([index_col_name, "sample_path", "group_name"])
        digest_keeper.single_genes_ds = pd.concat([digest_keeper.single_genes_ds, single_genes_total_ds], axis=0)
        gene_groups_ds = gene_groups_total_df.reset_index().melt(id_vars=["keyword"], var_name="sample_path").rename(columns={"value": digest_value_col_name})
        gene_groups_ds["group_name"] = group_name
        gene_groups_ds = gene_groups_ds.set_index(["keyword", "sample_path", "group_name"])
        digest_keeper.gene_groups_ds = pd.concat([digest_keeper.gene_groups_ds, gene_groups_ds], axis=0)
    digest_keepers_dict[digest_value_col_name] = digest_keeper

# Prepare p-values count
raw_ds = pd.DataFrame()
processed_ds = pd.DataFrame()
for digest_value_col_name in digest_value_col_names:
    digest_keeper = digest_keepers_dict.get(digest_value_col_name)
    raw_ds = pd.concat([raw_ds, digest_keeper.single_genes_ds], axis=1)
    processed_ds = pd.concat([processed_ds, digest_keeper.gene_groups_ds], axis=1)
    groupdata_file = "{}{}.groupdata".format(digest_keeper.samples_dir, "_".join(group_names))
    digest_keeper.groupdata_df.to_csv(groupdata_file, sep="\t", index=None, header=None)
    guideliner = GroupDataAssemblyGuideLiner(groupdata=groupdata_file,
                                             index_column="keyword",
                                             value_column=digest_value_col_name,
                                             output_dir="{}pvals/".format(digest_keeper.samples_dir))
    print("{}\n".format(guideliner.external_launch_command))

"""
# Pre-setup to launch from different node for group data file "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/two_first_checkpoints/gene_groups/log2(RPM+1)/1st_2nd_control.groupdata" and value column "log2(RPM+1)":

export IMG=ivasilyev/curated_projects:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it $IMG python3 /home/docker/scripts/statistical_tools/groupdata2statistics.py \
-g "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/two_first_checkpoints/gene_groups/log2(RPM+1)/1st_2nd_control.groupdata" \
-p "" \
-s "" \
-i "keyword" \
-v "log2(RPM+1)" \
-o "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/two_first_checkpoints/gene_groups/log2(RPM+1)/pvals/"



# Pre-setup to launch from different node for group data file "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/two_first_checkpoints/gene_groups/kRPKM/1st_2nd_control.groupdata" and value column "kRPKM":

export IMG=ivasilyev/curated_projects:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it $IMG python3 /home/docker/scripts/statistical_tools/groupdata2statistics.py \
-g "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/two_first_checkpoints/gene_groups/kRPKM/1st_2nd_control.groupdata" \
-p "" \
-s "" \
-i "keyword" \
-v "kRPKM" \
-o "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/two_first_checkpoints/gene_groups/kRPKM/pvals/"
"""

# Visualize boxplot data
for digest_value_col_name in digest_value_col_names:
    sns.set(style="whitegrid", font_scale=0.5)
    multiboxplot_dir = "{a}gene_groups/{b}/multiboxplots/".format(a=outputDir, b=digest_value_col_name)
    os.makedirs(multiboxplot_dir, exist_ok=True)
    fig, axes = plt.subplots(nrows=3, ncols=6, figsize=(10, 5), sharey=False)
    multiboxplot_df = processed_ds.reset_index().set_index("keyword")
    for ax, keyword in zip(axes.flatten(), sorted(DigestAssociationsKeeper.DRUG_CLASSES) + sorted(DigestAssociationsKeeper.RESISTANCE_MECHANISMS)):
        multiboxplot_data = multiboxplot_df.loc[keyword, :]
        multiboxplot_data.reset_index().rename(columns={keyword: digest_value_col_name}).to_csv("{a}data_{b}_{c}.tsv".format(a=multiboxplot_dir, b=digest_value_col_name, c=keyword), sep="\t", header=True, index=False)
        sns.boxplot(x="keyword", y=digest_value_col_name, hue="group_name", data=multiboxplot_data.reset_index(), orient="v", fliersize=1, linewidth=1, palette="Set3", ax=ax)
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc="right", bbox_to_anchor=(0.975, 0.5), title="Group ID", fancybox=True)
        ax.legend_.remove()
        ax.set_title(keyword.replace(" ", "\n"))
        ax.title.set_position([0.5, 0.97])
        ax.axes.get_xaxis().set_visible(False)
        ax.yaxis.label.set_visible(False)
        ax.tick_params(axis="y", which="major", labelrotation=0, pad=-3)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    ax0 = fig.add_axes([0, 0, 1, 1])
    plt.text(0.09, 0.5, digest_value_col_name, horizontalalignment="left", verticalalignment='center', rotation=90,
             transform=ax0.transAxes)
    plt.text(0.5, 0.95, "The abundance of ARGs based on time checkpoint groups after H.pylori eradication therapy",
             horizontalalignment="center", verticalalignment='center', transform=ax0.transAxes,
             fontsize="large", fontstyle="normal", fontweight="bold")
    ax0.set_axis_off()
    multiboxplot_image = "{a}multiboxplot_{b}.png".format(a=multiboxplot_dir, b=digest_value_col_name)
    fig.savefig(multiboxplot_image, format="png", dpi=900)
    plt.clf()
    plt.close()
