#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3
"""

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

# Prepare new raw reads for filtering to HG19 DB
projectDescriber = ProjectDescriber()
projectDescriber.sampledata = "/data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints_no_hg19.sampledata"
raw_reads_dirs = ["/data1/bio/170922/",
                  "/data1/bio/180419_NB501097_0017_AHL55TBGX5/HP/",
                  "/data1/bio/181017_NB501097_0024_AHL553BGX5/Conversion/"]


def create_sampledata_dict(dirs: list):
    output_dict = {}
    for raw_reads_dir in dirs:
        raw_reads_dir = Utilities.ends_with_slash(raw_reads_dir)
        files_list = os.listdir(raw_reads_dir)
        for file_name in files_list:
            if any([file_name.endswith(i) for i in ["csfasta", "fasta", "fa", "fastq", "fq", "gz"]]):
                sample_name = file_name.split("_")[0].strip()
                file_extension = file_name.split(".")[-1]
                sample_files = sorted(sorted([raw_reads_dir + i for i in files_list if len(re.findall("^{}".format(sample_name), i)) > 0 and i.endswith(file_extension)], key=len, reverse=True)[:2])
                sample_name = re.sub("_+", "_", re.sub("[^A-Za-z0-9]+", "_", sample_name))
                existing_files = output_dict.get(sample_name)
                if not existing_files:
                    output_dict[sample_name] = sample_files
                else:
                    if existing_files[0].endswith("csfasta") and file_extension != "csfasta":
                        print("Replacing SOLID reads with Illumina reads for sample {a}: {b}, {c}".format(a=sample_name, b=existing_files, c=sample_files))
                        output_dict[sample_name] = sample_files
    output_dict = dict(sorted(output_dict.items()))
    return output_dict


# Create sampledata for Illumina raw reads
raw_reads_dict = create_sampledata_dict(raw_reads_dirs)
raw_reads_dict = {k: raw_reads_dict[k] for k in raw_reads_dict if "HP" in k and not os.path.isfile("/data2/bio/Metagenomes/HG19/Unmapped_reads/{}_no_hg19.1.gz".format(k))}

Utilities.dump_2d_array([[k] + raw_reads_dict[k] for k in raw_reads_dict], file=projectDescriber.sampledata)

# Prepare deploy charts
launchGuideLiner = LaunchGuideLiner(charts_dir="{}{}/charts/".format(Utilities.ends_with_slash(projectDescriber.directory), "hg19"),
                                    deploy_prefix=projectDescriber.owner + "-hg19",
                                    nodes_number=7,
                                    threads_number="half",
                                    sampledata_file=projectDescriber.sampledata,
                                    refdata_file="/data/reference/homo_sapiens/hg/hg19/hg19.refdata",
                                    output_mask="hg19",
                                    output_dir="/data2/bio/Metagenomes/HG19")
launchGuideLiner.generate_config()
launchGuideLiner.get_deploy_guide()

"""
# Charts directory: '/data1/bio/projects/dsafina/hp_checkpoints/hg19/charts/'

# Look for Redis pod & service:
kubectl get pods

# Deploy master chart
kubectl create -f /data1/bio/projects/dsafina/hp_checkpoints/hg19/charts/master.yaml

# Deploy worker chart only when the master pod ('dsafina-hg19-queue') is complete
kubectl create -f /data1/bio/projects/dsafina/hp_checkpoints/hg19/charts/worker.yaml

# View active nodes
kubectl describe pod dsafina-hg19-job- | grep Node:

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod dsafina-hg19-queue
kubectl delete job dsafina-hg19-job
"""

# Disable coverage extraction
subprocess.getoutput("sed -i 's|/home/docker/scripts/master.py, -s, |/home/docker/scripts/master.py, -n, -s, |g' /data1/bio/projects/dsafina/hp_checkpoints/hg19/charts/master.yaml")

# Deploy from MASTER

# Prepare sample data for CARD DB
filtered_reads_dirs = ["/data2/bio/Metagenomes/HG19/Non-mapped_reads/",
                       "/data2/bio/Metagenomes/HG19/Unmapped_reads/"]

filtered_reads_dict = create_sampledata_dict(filtered_reads_dirs)
filtered_reads_dict = {k: filtered_reads_dict[k] for k in filtered_reads_dict if any(i in k for i in ["HP", "I", "SRR"])}
Utilities.dump_2d_array([[k] + raw_reads_dict[k] for k in raw_reads_dict], file=projectDescriber.sampledata)

"""
# Remove already mapped sample data

export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/verify_coverages.py \
-i /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints_no_hg19.sampledata \
-p /data2/bio/Metagenomes/CARD/Statistics/ \
-s _card_v2.0.3_coverage.tsv \
-g /data/reference/CARD/card_v2.0.3/index/card_v2.0.3_samtools.genome \
-d
"""

"""
Done.
Files to process: 17
Dumped sample data: '/data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-10-25-21-01-57.sampledata'
Dumped debug table: '/data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-10-25-21-01-57.sampledata_debug.tsv'
"""

referenceDescriber = ReferenceDescriber()

launchGuideLiner = LaunchGuideLiner(charts_dir="{}{}/charts/".format(Utilities.ends_with_slash(projectDescriber.directory), referenceDescriber.alias),
                                    deploy_prefix="{a}-{b}-{c}".format(a=projectDescriber.owner, b=projectDescriber.name, c=referenceDescriber.alias),
                                    nodes_number=7,
                                    threads_number="half",
                                    sampledata_file="/data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-10-25-21-01-57.sampledata",
                                    refdata_file=referenceDescriber.refdata,
                                    output_mask=referenceDescriber.alias,
                                    output_dir="/data2/bio/Metagenomes/{}/".format(referenceDescriber.name))
launchGuideLiner.generate_config()
launchGuideLiner.get_deploy_guide()

"""
# Charts directory: '/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/charts/'

# Look for Redis pod & service:
kubectl get pods

# (Optional) If required, flush all Redis queues:
export IMG=ivasilyev/bwt_filtering_pipeline_master:latest && \
docker pull $IMG && docker run --rm -it $IMG redis-cli -h redis flushall

# (Optional) Or remove Redis server by itself:
kubectl delete svc,pod redis

# (Optional) Deploy if not present:
kubectl create -f  https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/meta/charts/redis.yaml

# Deploy master chart
kubectl create -f /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/charts/master.yaml

# Deploy worker chart only when the master pod ('dsafina-hp-checkpoints-card-v2-0-3-queue') is complete
kubectl create -f /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/charts/worker.yaml

# View active nodes
kubectl describe pod dsafina-hp-checkpoints-card-v2-0-3-job- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/CARD/Statistics/*coverage.txt | wc -l) OF $(cat /data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-10-25-21-01-57.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod dsafina-hp-checkpoints-card-v2-0-3-queue
kubectl delete job dsafina-hp-checkpoints-card-v2-0-3-job

# Checkout (from WORKER node)
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/verify_coverages.py \
-i /data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-10-25-21-01-57.sampledata \
-p /data2/bio/Metagenomes/CARD/Statistics/ \
-s _card_v2.0.3_coverage.tsv \
-g /data/reference/CARD/card_v2.0.3/index/card_v2.0.3_samtools.genome \
-d
"""

"""
Done.
Files to process: 3
Dumped sample data: '/data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-10-26-06-41-34.sampledata'
Dumped debug table: '/data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-10-26-06-41-34.sampledata_debug.tsv'
"""

"""
# Again, manual launch for the small sample data:

export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/nBee.py \
-i /data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-10-26-06-41-34.sampledata \
-r /data/reference/CARD/card_v2.0.3/index/card_v2.0.3_refdata.json \
-m card_v2.0.3 -t half -o /data2/bio/Metagenomes/CARD/

# Checkout (from WORKER node)
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/verify_coverages.py \
-i /data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-10-26-06-41-34.sampledata \
-p /data2/bio/Metagenomes/CARD/Statistics/ \
-s _card_v2.0.3_coverage.tsv \
-g /data/reference/CARD/card_v2.0.3/index/card_v2.0.3_samtools.genome \
-d
"""

"""
Done.
Files to process: 0
"""

# Create new group data for 3 paired checkpoints
time_points_df = pd.read_excel("/data2/bio/Chocophlan/HP_samples.xlsx", sheet_name="3 time points", header=0)
time_points_df.dropna(thresh=2, inplace=True)
time_points_df.dropna(axis="columns", inplace=True)
time_points_df = time_points_df[list(time_points_df)].astype(int)
time_points_df = time_points_df.melt(var_name="group_name", value_name="sample_name")
time_points_df["sample_name"] = time_points_df["sample_name"].astype(str) + "HP"

time_points_df.loc[:, ["sample_name", "group_name"]].to_csv("/data1/bio/projects/dsafina/hp_checkpoints/3_paired_hp_checkpoints.groupdata", sep='\t', header=False, index=False)

# os.rename("/data2/bio/Metagenomes/CARD/Statistics/216_1HP_card_v2.0.3_coverage.tsv", "/data2/bio/Metagenomes/CARD/Statistics/216HP_card_v2.0.3_coverage.tsv")

"""
# Combine data for RPM
rm -rf /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/pvals /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/metadata_digest

export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git clone https://github.com/ivasilyev/statistical_tools.git
cd statistical_tools

python3 groupdata2statistics.py -g /data1/bio/projects/dsafina/hp_checkpoints/3_paired_hp_checkpoints.groupdata \
-p /data2/bio/Metagenomes/CARD/Statistics/ \
-s _card_v2.0.3_coverage.tsv \
-i reference_id \
-v id_mapped_reads_per_million_sample_mapped_reads \
-o /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/pvals/RPM/
"""

"""
# Combine data for RPKM

export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git clone https://github.com/ivasilyev/statistical_tools.git
cd statistical_tools

python3 groupdata2statistics.py -g /data1/bio/projects/dsafina/hp_checkpoints/3_paired_hp_checkpoints.groupdata \
-p /data2/bio/Metagenomes/CARD/Statistics/ \
-s _card_v2.0.3_coverage.tsv \
-i reference_id \
-v id_mapped_reads_per_kbp_per_million_sample_mapped_reads \
-o /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/pvals/RPKM/
"""

index_col_name = "reference_id"
annotation_df = pd.read_table("/data/reference/CARD/card_v2.0.3/index/card_v2.0.3_annotation.tsv").set_index(index_col_name)


class ReversedGroupComparator:
    def __init__(self, name):
        self.name = str(name)
        self.columns = []
        self.rows = []
        self.comparison_2d_array = []
        self.comparisons = []


drug_classes_dict = {"beta-lactam": ("cephalosporin", "penam", "penem"),
                     "aminoglycoside": (),
                     "fluoroquinolone": (),
                     "glycopeptide antibiotic": (),
                     "lincosamide": (),
                     "macrolide": (),
                     "nucleoside antibiotic": (),
                     "peptide antibiotic": (),
                     "phenicol": (),
                     "sulfonamide": (),
                     "tetracycline": (),
                     "triclosan": ()}
resistance_mechanisms_dict = {"efflux",
                              "inactivation",
                              "reduced permeability",
                              "target alteration",
                              "target protection",
                              "target replacement"}

keywords_list = sorted(list(drug_classes_dict)) + sorted(list(resistance_mechanisms_dict))
group_names = set(time_points_df["group_name"])

summarized_dss_dict = {}
for pivot_value_col_name in ("RPM", "RPKM"):
    total_df_file = "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/pvals/{}/1st_2nd_3rd_total_dataframe.tsv".format(pivot_value_col_name)
    total_df = pd.read_table(total_df_file).set_index(index_col_name)
    annotated_total_df = pd.concat([annotation_df, total_df], axis=1)
    annotated_total_df.index.name = index_col_name
    annotated_total_df.to_csv(total_df_file.replace("_total_dataframe.tsv", "_total_dataframe_annotated.tsv"), sep='\t', header=True, index=True)
    #
    groupdata_comparators_dict = {}
    for first_group_name in group_names:
        reversedGroupComparator = ReversedGroupComparator(first_group_name)
        for annotated_total_df_column_name in list(annotated_total_df):
            annotated_total_df_column_name_digest_list = annotated_total_df_column_name.split("_")
            if annotated_total_df_column_name_digest_list[0].strip() == first_group_name and len(annotated_total_df_column_name_digest_list) > 0 and annotated_total_df_column_name_digest_list[1].strip().startswith("/"):
                reversedGroupComparator.columns.append(annotated_total_df_column_name)
        reversedGroupComparator.rows = annotated_total_df.loc[annotated_total_df["{}_non-zero_values".format(first_group_name)].astype(int) > 0, :].index.tolist()
        for second_group_name in [i for i in group_names if i != first_group_name]:
            comparison_pairs = ((first_group_name, second_group_name), (second_group_name, first_group_name))
            for comparison_pair in comparison_pairs:
                comparison_name = "{a}_vs_{b}_is_rejected_by_fdr_bh_for_wilcoxon".format(a=comparison_pair[0], b=comparison_pair[1])
                if comparison_name in list(annotated_total_df):
                    reversedGroupComparator.comparisons.append(comparison_name)
                    reversedGroupComparator.comparison_2d_array.append(comparison_pair)
        groupdata_comparators_dict[first_group_name] = reversedGroupComparator
    #
    values_columns_list = [j for i in [groupdata_comparators_dict[k].columns for k in groupdata_comparators_dict] for j in i]
    #
    boxplot_dfs_dict = {}
    for drug_class in drug_classes_dict:
        if len(drug_classes_dict[drug_class]) > 0:
            boxplot_dfs_dict[drug_class] = annotated_total_df.loc[annotated_total_df["Drug Class"].apply(lambda x: any(i in str(x) for i in drug_classes_dict[drug_class])) == True, values_columns_list]
        else:
            boxplot_dfs_dict[drug_class] = annotated_total_df.loc[annotated_total_df["Drug Class"].str.contains(drug_class) == True, values_columns_list]
    #
    for resistance_mechanism in resistance_mechanisms_dict:
        boxplot_dfs_dict[resistance_mechanism] = annotated_total_df.loc[annotated_total_df["Resistance Mechanism"].str.contains(resistance_mechanism) == True, values_columns_list]
    #
    summarized_series_list = [boxplot_dfs_dict[k].sum().rename(k) for k in boxplot_dfs_dict]
    #
    summarized_df = pd.concat(summarized_series_list, axis=1)
    summarized_df.index.name = "coverage_file"
    summarized_df = summarized_df.transpose().fillna(0)
    summarized_df.index.name = "keywords"
    #
    reversed_groupdata_2d_array = []
    for col_name in list(summarized_df):
        group_name = col_name.split("_")[0]
        coverage_file = "_".join(col_name.split("_")[1:])
        sample_name = re.sub("[\W]+", "_", coverage_file.replace("/data2/bio/Metagenomes/CARD/Statistics/", "").replace("_card_v2.0.3_coverage.tsv", ""))
        isolated_value_dir = "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/metadata_digest/{a}/{b}/".format(a=pivot_value_col_name, b=group_name)
        os.makedirs(isolated_value_dir, exist_ok=True)
        isolated_value_file = "{}{}.tsv".format(isolated_value_dir, sample_name)
        summarized_df[col_name].reset_index().rename(columns={col_name: pivot_value_col_name}).to_csv(isolated_value_file, sep='\t', header=True, index=False)
        reversed_groupdata_2d_array.append([isolated_value_file, group_name])
    #
    Utilities.dump_2d_array(array=reversed_groupdata_2d_array, file="/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/metadata_digest/{a}/{a}.groupdata".format(a=pivot_value_col_name))
    #
    summarized_ds = pd.melt(summarized_df.reset_index(), id_vars=["keywords"], var_name="coverage_file", value_name=pivot_value_col_name).set_index(["keywords", "coverage_file"])
    summarized_dss_dict[pivot_value_col_name] = summarized_ds

digest_values_ds = pd.concat([summarized_dss_dict[k] for k in summarized_dss_dict], axis=1).reset_index()

digest_values_ds["group_name"] = digest_values_ds["coverage_file"].apply(lambda x: x.split("_")[0])
digest_values_ds["log2(RPM+1)"] = np.log2(digest_values_ds["RPM"] + 1)
digest_values_ds["kRPKM"] = digest_values_ds["RPKM"].astype(float) / 1000.0

digest_values_ds.to_csv("/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/metadata_digest/digest_values_dataset.tsv", sep="\t", header=True, index=False)

# digest_values_ds = pd.read_table("/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/metadata_digest/digest_values_dataset.tsv", sep="\t", header=0, engine="python")

"""
# Combine digested data for RPM

export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git clone https://github.com/ivasilyev/statistical_tools.git
cd statistical_tools

python3 groupdata2statistics.py -g /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/metadata_digest/RPM/RPM.groupdata \
-i keywords \
-v RPM \
-o /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/metadata_digest/pvals/RPM
"""

"""
# Combine digested data for RPKM

export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git clone https://github.com/ivasilyev/statistical_tools.git
cd statistical_tools

python3 groupdata2statistics.py -g /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/metadata_digest/RPKM/RPKM.groupdata \
-i keywords \
-v RPKM \
-o /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/metadata_digest/pvals/RPKM
"""

# Visualize boxplot data
for boxplot_y_col_name in ("log2(RPM+1)", "kRPKM"):
    sns.set(style="whitegrid", font_scale=0.5)
    multiboxplot_alias = re.sub("[\W]+", "_", boxplot_y_col_name).strip("_")
    multiboxplot_dir = "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/metadata_digest/multiboxplots/{}/".format(multiboxplot_alias)
    os.makedirs(os.path.dirname(multiboxplot_dir), exist_ok=True)
    #
    fig, axes = plt.subplots(nrows=3, ncols=6, figsize=(10, 5), sharey=False)
    for ax, keyword in zip(axes.flatten(), keywords_list):
        multiboxplot_data = digest_values_ds.loc[digest_values_ds["keywords"] == keyword, ["keywords", boxplot_y_col_name, "group_name"]]
        multiboxplot_data.to_csv("{a}dataset_{b}_{c}.tsv".format(a=multiboxplot_dir, b=multiboxplot_alias, c=keyword), sep="\t", header=True, index=False)
        sns.boxplot(x="keywords", y=boxplot_y_col_name, hue="group_name", data=multiboxplot_data, orient="v", fliersize=1, linewidth=1, palette="Set3", ax=ax)
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
    plt.text(0.09, 0.5, boxplot_y_col_name, horizontalalignment="left", verticalalignment='center', rotation=90,
             transform=ax0.transAxes)
    plt.text(0.5, 0.95, "The abundance of ARGs based on time checkpoint groups after H.pylori eradication therapy",
             horizontalalignment="center", verticalalignment='center', transform=ax0.transAxes,
             fontsize="large", fontstyle="normal", fontweight="bold")
    ax0.set_axis_off()
    multiboxplot_image = "{a}multiboxplot_{b}.png".format(a=multiboxplot_dir, b=multiboxplot_alias)
    fig.savefig(multiboxplot_image, format="png", dpi=900)
    plt.clf()
    plt.close()
