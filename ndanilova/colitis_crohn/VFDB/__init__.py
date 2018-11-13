#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -e DISPLAY=$DISPLAY -it ${DOCKER_IMAGE_NAME} python3
"""

from ndanilova.colitis_crohn.ProjectDescriber import ProjectDescriber
from meta.scripts.vfdb.ReferenceDescriber import ReferenceDescriber
from meta.scripts.LaunchGuideLiner import LaunchGuideLiner
from meta.scripts.Utilities import Utilities
from meta.scripts.GroupDataPreparer import GroupDataAssemblyGuideLiner
import subprocess
import os
import re
import pandas as pd
import json
from collections import Counter


from meta.scripts.PivotTableAnnotator import PivotTableAnnotator


outputDir = "/data1/bio/projects/ndanilova/colitis_crohn/VFDB/"
os.makedirs(outputDir, exist_ok=True)

guideLiner = LaunchGuideLiner(charts_dir="/data1/bio/projects/ndanilova/colitis_crohn/VFDB/charts/",
                              deploy_prefix="ndanilova-vfdb",
                              nodes_number=9,
                              threads_number="half",
                              sampledata_file=ProjectDescriber.sampledata,
                              refdata_file=ReferenceDescriber.refdata,
                              output_mask="no_hg19",
                              output_dir="/data2/bio/Metagenomes/Toxins/VFDB/")

guideLiner.generate_config()
guideLiner.get_deploy_guide()

"""
# Charts directory: '/data1/bio/projects/ndanilova/colitis_crohn/VFDB/charts/'

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
kubectl create -f /data1/bio/projects/ndanilova/colitis_crohn/VFDB/charts/master.yaml

# Deploy worker chart only when the master pod ('ndanilova-vfdb-queue') is complete
kubectl create -f /data1/bio/projects/ndanilova/colitis_crohn/VFDB/charts/worker.yaml

# View active nodes
kubectl describe pod ndanilova-vfdb-job- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/Toxins/VFDB/Statistics/*coverage.txt | wc -l) OF $(cat /data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod ndanilova-vfdb-queue
kubectl delete job ndanilova-vfdb-job

# Checkout (from WORKER node)
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/verify_coverages.py \
-i /data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.sampledata \
-p /data2/bio/Metagenomes/Toxins/VFDB/Statistics/ \
-s _vfdb_v2018.11.09_coverage.tsv \
-g /data/reference/VFDB/vfdb_v2018.11.09/index/vfdb_v2018.11.09_samtools.genome \
-d
"""

"""
Files to process: 0
Dumped sample data: /data2/bio/Metagenomes/Toxins/VFDB/Statistics/sampledata/2018-11-10-12-27-21.sampledata
Dumped debug table: '/data2/bio/Metagenomes/Toxins/VFDB/Statistics/sampledata/2018-11-10-12-27-21.sampledata_debug.tsv'
"""

index_col_name = "reference_id"
value_col_names_abbreviations_dict = {"id_mapped_reads_per_million_sample_mapped_reads": "RPM"}
prefix = "/data2/bio/Metagenomes/Toxins/VFDB/Statistics/"
suffix = "_vfdb_v2018.11.09_coverage.tsv"

# TODO: Add RPKM

class DanilovaAwesomeGroupAnalysisHandler:
    def __init__(self, index_col_name, value_col_name, value_col_name_abbreviation):
        # Fields to fill before first data assembly for raw coverage files with gene name as index
        self.index_col_name = index_col_name
        self.pivot_value_col_name = value_col_name
        self.pivot_value_col_name_abbreviation = value_col_name_abbreviation
        self.groupdata_file = ""
        self.groupdata_digest_name = ""
        self.groupdata_dict = {}
        self.raw_all_sample_names_list = []
        self.output_dir = ""
        self.raw_pvals_dir = ""
        self.raw_pvals_guideline = ""
        self.raw_prefix = ""
        self.raw_suffix = ""
        # Fields to fill before the second data assembly:
        self.raw_pvals_df = pd.DataFrame()
        self.raw_pivot_dfs_dict = {}
        self.annotation_df = pd.DataFrame()
        self.digest_virulence_dfs_dict = {}
        # Fields to collect data and start the second data assembly:
        self.digest_dir = ""
        # Virulence:
        self.digest_virulence_pvals_dir = ""
        self.digest_virulence_samples_dir = ""
        self.digest_virulence_groupdata = ""
        self.digest_virulence_guidelines = ""
        # Genera:
        self.digest_genera_groupdata = ""
        self.digest_genera_guidelines = ""
        # Fields required for visualization
        self.digest_virulence_pivot_df = pd.DataFrame()
    def set_groupdata_dict(self, groupdata_file: str):
        self.groupdata_file = groupdata_file
        self.groupdata_digest_name = Utilities.filename_only(self.groupdata_file).replace(".groupdata", "")
        groupdata_df = pd.read_table(self.groupdata_file, sep="\t", header="infer", names=["sample_name", "group_name"])
        self.groupdata_dict = {i: sorted(set(
            groupdata_df.loc[groupdata_collection_df["group_name"] == i, ["sample_name"]])) for i in sorted(
            set(groupdata_df["group_name"]))}
        self.raw_all_sample_names_list = sorted(set(groupdata_df["sample_name"]))
    def set_raw_pvals_dir(self, output_dir: str):
        self.output_dir = Utilities.ends_with_slash(output_dir)
        self.raw_pvals_dir = "{OUTPUT_DIR}pvals/{VALUE_COLUMN}/{GROUPS}/".format(OUTPUT_DIR=self.output_dir,
                                                                                 VALUE_COLUMN=self.pivot_value_col_name_abbreviation,
                                                                                 GROUPS=self.groupdata_digest_name)
        self.digest_dir = "{OUTPUT_DIR}digest/{VALUE_COLUMN}/{GROUPS}/".format(OUTPUT_DIR=self.output_dir,
                                                                               VALUE_COLUMN=self.pivot_value_col_name_abbreviation,
                                                                               GROUPS=self.groupdata_digest_name)
    def get_raw_guideliners(self, prefix, suffix):
        self.raw_prefix = prefix
        self.raw_suffix = suffix
        guideliner = GroupDataAssemblyGuideLiner(groupdata=self.groupdata_file,
                                                 prefix=self.raw_prefix,
                                                 suffix=self.raw_suffix,
                                                 index_column=self.index_col_name,
                                                 value_column=self.pivot_value_col_name,
                                                 output_dir=self.raw_pvals_dir)
        self.raw_pvals_guideline = guideliner.external_launch_command
        print(self.raw_pvals_guideline)
    def _load_df(self, file: str):
        return pd.read_table(file, sep="\t", header=0).set_index(self.index_col_name).sort_index()
    def set_assembled_raw_data_paths(self):
        completed_raw_pvals_file = subprocess.getoutput("ls -d {}*_total_dataframe.tsv".format(self.raw_pvals_dir)).strip()
        self.raw_pvals_df = self._load_df(completed_raw_pvals_file)
        # self.raw_pivot_files_list = [i.strip() for i in subprocess.getoutput("ls -d {}pivot_by_*".format(self.raw_pvals_dir)).split("\n")]
        self.raw_pivot_dfs_dict = {k: self._load_df("{PVALS_DIR}pivot_by_{VALUE_COLUMN}_{GROUP}.tsv".format(
            PVALS_DIR=self.raw_pvals_dir, VALUE_COLUMN=self.pivot_value_col_name, GROUP=k)) for k in self.groupdata_dict}
    def annotate_assembled_raw_data(self, annotation_file: str):
        self.annotation_df = self._load_df(annotation_file)
        # Annotate pvals table
        self.raw_pvals_df = pd.concat([self.annotation_df, self.raw_pvals_df], axis=1, sort=True)
        self.raw_pvals_df.index.names = [self.index_col_name]
        # Annotate pivots
        for group in self.raw_pivot_dfs_dict:
            df = self.raw_pivot_dfs_dict[group]
            df = pd.concat([self.annotation_df, df], axis=1, sort=True)
            df.index.names = [self.index_col_name]
            self.raw_pivot_dfs_dict[group] = df
    def digest_virulence_pivots(self):
        groupdata_2d_array = []
        self.digest_virulence_samples_dir = "{}virulence/samples/".format(self.digest_dir)
        search_col_names = ["gene_description",]
        for group in self.raw_pivot_dfs_dict:
            raw_df = self.raw_pivot_dfs_dict[group].loc[:,
                     search_col_names + [i for i in list(self.raw_pivot_dfs_dict[group]) if i not in list(self.annotation_df)]]
            digest_df = DigestAssociationsKeeper.digest_df(raw_df, DigestAssociationsKeeper.VIRULENCE_FACTORS, search_col_names)
            self.digest_virulence_dfs_dict[group] = digest_df
            for sample_path in list(digest_df):
                sample_name = sample_path.replace(self.raw_prefix, "").replace(self.raw_suffix, "")
                self.digest_virulence_samples_dir = "{}virulence/samples/{}/".format(self.digest_dir, group)
                os.makedirs(self.digest_virulence_samples_dir, exist_ok=True)
                output_file = "{}{}.tsv".format(self.digest_virulence_samples_dir, sample_name)
                digest_df[sample_path].reset_index().rename(
                    columns={sample_path: self.pivot_value_col_name}).to_csv(output_file, sep='\t', header=True, index=False)
                groupdata_2d_array[group].append([group, output_file])
        self.digest_virulence_pvals_dir = "{}virulence/pvals/".format(self.digest_dir)
        os.makedirs(self.digest_virulence_pvals_dir, exist_ok=True)
        self.digest_virulence_groupdata = "{}{}.groupdata".format(self.digest_virulence_pvals_dir, "_".join(list(self.raw_pivot_dfs_dict)))
        Utilities.dump_2d_array(groupdata_2d_array, file=self.digest_virulence_groupdata)
    def get_digest_virulence_guideliners(self):
        guideliner = GroupDataAssemblyGuideLiner(groupdata=self.digest_virulence_groupdata,
                                                 prefix=self.raw_prefix,
                                                 suffix=self.raw_suffix,
                                                 index_column=self.index_col_name,
                                                 value_column=self.pivot_value_col_name,
                                                 output_dir=self.digest_virulence_pvals_dir)
        self.raw_pvals_guideline = guideliner.external_launch_command
        print(self.raw_pvals_guideline)


# First run
handlers_dict = {}
for pivot_value_col_name in value_col_names_abbreviations_dict:
    for groupdata_file in ProjectDescriber.groupdata:
        handler = DanilovaAwesomeGroupAnalysisHandler(index_col_name="reference_id",
                                                      value_col_name="id_mapped_reads_per_million_sample_mapped_reads",
                                                      value_col_name_abbreviation="RPM")
        handler.set_groupdata_dict(groupdata_file)
        handler.set_raw_pvals_dir(outputDir)
        handler.get_raw_guideliners(prefix, suffix)
        # Note the reversed keys order
        if not handlers_dict[handler.groupdata_digest_name]:
            handlers_dict[handler.groupdata_digest_name] = {}
        handlers_dict[handler.groupdata_digest_name][pivot_value_col_name] = handler

# Data assembly for per-gene p-values count
"""
# Pre-setup to launch from different node for group data file "/data1/bio/projects/ndanilova/colitis_crohn/group_data_digest/all_groups.groupdata" and value column "id_mapped_reads_per_million_sample_mapped_reads":

export IMG=ivasilyev/curated_projects:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it $IMG python3 /home/docker/statistical_tools/groupdata2statistics.py \
-g /data1/bio/projects/ndanilova/colitis_crohn/group_data_digest/all_groups.groupdata \
-p /data2/bio/Metagenomes/Toxins/VFDB/Statistics/ \
-s _vfdb_v2018.11.09_coverage.tsv \
-i reference_id \
-v id_mapped_reads_per_million_sample_mapped_reads \
-o /data1/bio/projects/ndanilova/colitis_crohn/VFDB/pvals/RPM/all_groups/


# Pre-setup to launch from different node for group data file "/data1/bio/projects/ndanilova/colitis_crohn/group_data_digest/remission_vs_escalation.groupdata" and value column "id_mapped_reads_per_million_sample_mapped_reads":

export IMG=ivasilyev/curated_projects:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it $IMG python3 /home/docker/statistical_tools/groupdata2statistics.py \
-g /data1/bio/projects/ndanilova/colitis_crohn/group_data_digest/remission_vs_escalation.groupdata \
-p /data2/bio/Metagenomes/Toxins/VFDB/Statistics/ \
-s _vfdb_v2018.11.09_coverage.tsv \
-i reference_id \
-v id_mapped_reads_per_million_sample_mapped_reads \
-o /data1/bio/projects/ndanilova/colitis_crohn/VFDB/pvals/RPM/remission_vs_escalation/


# Pre-setup to launch from different node for group data file "/data1/bio/projects/ndanilova/colitis_crohn/group_data_digest/crohn_vs_colitis.groupdata" and value column "id_mapped_reads_per_million_sample_mapped_reads":

export IMG=ivasilyev/curated_projects:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it $IMG python3 /home/docker/statistical_tools/groupdata2statistics.py \
-g /data1/bio/projects/ndanilova/colitis_crohn/group_data_digest/crohn_vs_colitis.groupdata \
-p /data2/bio/Metagenomes/Toxins/VFDB/Statistics/ \
-s _vfdb_v2018.11.09_coverage.tsv \
-i reference_id \
-v id_mapped_reads_per_million_sample_mapped_reads \
-o /data1/bio/projects/ndanilova/colitis_crohn/VFDB/pvals/RPM/crohn_vs_colitis/
"""

# Second run
for groupdata_digest in handlers_dict:
    for pivot_value_col_name in handlers_dict[groupdata_digest]:
        handler = handlers_dict[groupdata_digest][pivot_value_col_name]
        handler.set_assembled_raw_data_paths()
        handler.annotate_assembled_raw_data("/data/reference/VFDB/vfdb_v2018.11.09/index/vfdb_v2018.11.09_annotation.tsv")
        handler.digest_virulence_pivots()


# Get guidelines for data assembly for all group data
assembled_pvals_dirs_dict = {}
for pivot_value_col_name in value_col_names_abbreviations_dict:
    assembled_pvals_dirs_dict[pivot_value_col_name] = {}
    for groupdata_file in ProjectDescriber.groupdata:
        groupdata_digest_name = Utilities.filename_only(groupdata_file).replace(".groupdata", "")
        assembled_pvals_dir = "{OUTPUT_DIR}pvals/{PIVOT_VALUE_COLUMN}/{GROUPS}/".format(OUTPUT_DIR=outputDir,
                                                                                        PIVOT_VALUE_COLUMN=value_col_names_abbreviations_dict[pivot_value_col_name],
                                                                                        GROUPS=groupdata_digest_name)
        assembled_pvals_dirs_dict[pivot_value_col_name][groupdata_digest_name] = assembled_pvals_dir
        guideLiner = GroupDataAssemblyGuideLiner(groupdata=groupdata_file,
                                                 prefix=prefix,
                                                 suffix=suffix,
                                                 index_column=index_col_name,
                                                 value_column=pivot_value_col_name,
                                                 output_dir=assembled_pvals_dir)
        print(guideLiner.external_launch_command)

# Collect all group data
groupdata_collection_df = pd.concat([pd.read_table(i, header="infer", names=["sample_name", "group_name"]) for i in ProjectDescriber.groupdata for i in ProjectDescriber.groupdata], axis=0, ignore_index=True)
groupdata_digest_dict = {i: sorted(set(groupdata_collection_df.loc[groupdata_collection_df["group_name"] == "colitis", "sample_name"])) for i in sorted(set(groupdata_collection_df["group_name"]))}
sample_names = sorted(set(groupdata_collection_df["sample_name"]))

# Join and annotate it
annotation_df = pd.read_table("/data/reference/VFDB/vfdb_v2018.11.09/index/vfdb_v2018.11.09_annotation.tsv").set_index(index_col_name).sort_index()

assembled_pvals_dfs_dict = {}
for pivot_value_col_name in assembled_pvals_dirs_dict:
    assembled_pvals_dfs_dict[pivot_value_col_name] = {}
    group_digests_dict = assembled_pvals_dirs_dict[pivot_value_col_name]
    for group_digest in group_digests_dict:
        assembled_pvals_dir = group_digests_dict[group_digest]
        assembled_pvals_total_file = subprocess.getoutput("ls -d {}*_total_dataframe.tsv".format(assembled_pvals_dir))
        assembled_pvals_total_df = pd.concat([annotation_df.copy(),
                                              pd.read_table(assembled_pvals_total_file, sep="\t", header=0).set_index(
                                                  index_col_name).sort_index()], axis=1, sort=True)
        assembled_pvals_total_df.index.names = [index_col_name]
        assembled_pvals_dfs_dict[pivot_value_col_name][group_digest] = assembled_pvals_total_df


class DigestAssociationsKeeper:
    VIRULENCE_FACTORS = {"adhesion": ("adhesin", "adhesion", "laminin"),
                         "invasion": ("invasion", "invasin"),
                         "iron metabolism": (
                             "iron", "siderophore", "ferric", "ferrienterobactin", "ferrochelatase", "ferrichrome",
                             "aerobactin", "enterochelin", "enterobactin", "yersiniabactin", "yersinabactin",
                             "ferrienterobactin", "chrysobactin", "ornibactin", "precolibactin", "colibactin",
                             "ferripyoverdine"),
                         "pili-related": (
                             "pilin", "pilus", "prepilin", "fimbrial", "fimbriae", "fimbrillin", "fimbrin"),
                         "flagella-related": ("motor", "flagellin", "flagellar", "flagella", "flagellum"),
                         "regulation": ("regulator", "regulatory", "regulation", "receptor", "effector"),
                         "transport": ("transport", "transporter", "permease", "porin", "export"),
                         "toxin": (
                             "toxin", "enterotoxin", "cytotoxic", "necrotizing", "endotoxin", "leukotoxin", "exotoxin",
                             "cytolethal", "o-antigen", "lipooligosaccharide", "lipopolysaccharide"),
                         "chaperones": ("chaperone", "chaperonin"),
                         "cytolysins": ("hemolysin", "cytolysin"),
                         "cell wall related": ("cell wall",),
                         "cell membrane related": ("membrane",),
                         "capsule-related": ("capsular", "capsule"),
                         "lipoglycans": ("o-antigen", "lipooligosaccharide", "lipopolysaccharide"),
                         "secretion": ("secretion", "secreted", "secretory"),
                         "efflux": ("efflux",)}
    @staticmethod
    def digest_df(df: pd.DataFrame, associations: dict, columns_with_keywords: list):
        df_columns = list(df)
        columns_with_keywords = [i for i in columns_with_keywords if len(columns_with_keywords) > 0]
        if len(columns_with_keywords) == 0:
            raise ValueError("No column for keyword search specified")
        try:
            df["lookup_column"] = df.loc[:, columns_with_keywords].astype(str).apply(
                lambda x: " ".join([CounterWrapper.prepare_string(i) for i in x]), axis=1)
        except KeyError:
            print(list(df), associations, columns_with_keywords)
        keywords_series = []
        for main_word in associations:
            key_words = associations.get(main_word)
            if not key_words or len(key_words) == 0:
                key_words = (main_word, )
            sub_df = df.loc[df["lookup_column"].apply(lambda x: any(i.lower() in str(x) for i in key_words)) == True,
                            [i for i in df_columns if i not in columns_with_keywords]]
            keywords_series.append(sub_df.sum().rename(main_word))
        out_df = pd.concat(keywords_series, axis=1)
        out_df.index.name = "sample_path"
        out_df = out_df.transpose().fillna(0)
        out_df.index.name = "keyword"
        return out_df


class CounterWrapper:
    @staticmethod
    def prepare_string(s: str):
        return re.sub(" +", " ", re.sub("[^a-z0-9\-]+", " ", s.lower())).strip()
    @staticmethod
    def count_words_in_series(series: pd.Series):
        lst = [i for i in CounterWrapper.prepare_string(" ".join(series.fillna("").values.tolist())).split(" ") if len(i) > 3]
        return Counter(lst)
    @staticmethod
    def dump_counter(counter: Counter, file: str):
        Utilities.dump_2d_array([("keyword", "occurrences")] + counter.most_common(), file=file)


pivot_value_col_name = "RPM"
group_digest = "colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr"
total_df_file = "/data1/bio/projects/ndanilova/colitis_crohn/VFDB/pvals/{a}/{b}_total_dataframe.tsv".format(a=pivot_value_col_name, b=group_digest)
total_df = pd.read_table(total_df_file).set_index(index_col_name).sort_index()
annotated_total_df = pd.concat([annotation_df, total_df], axis=1, sort=True)
annotated_total_df.index.name = index_col_name
annotated_total_df.to_csv(total_df_file.replace("_total_dataframe.tsv", "_total_dataframe_annotated.tsv"), sep='\t',
                          header=True, index=True)

# Filter significant rows to get most common words
significant_pvals_df = annotated_total_df.loc[annotated_total_df["null_hypothesis_rejections_counter"].astype(int) > 0]

annotation_counter_col_name = "gene_description"

# Counting most common words in gene description
gene_description_counters = CounterWrapper.count_words_in_series(significant_pvals_df[annotation_counter_col_name])
CounterWrapper.dump_counter(gene_description_counters, file="/data1/bio/projects/ndanilova/colitis_crohn/VFDB/pvals/{a}/{b}_significant_{c}_counter.tsv".format(a=pivot_value_col_name, b=group_digest, c=annotation_counter_col_name))

# Static associations dict
gene_description_dict = {"adhesion": ("adhesin", "adhesion", "laminin"),
                         "invasion": ("invasion", "invasin"),
                         "iron metabolism": ("iron", "siderophore", "ferric", "ferrienterobactin", "ferrochelatase", "ferrichrome", "aerobactin", "enterochelin", "enterobactin", "yersiniabactin", "yersinabactin", "ferrienterobactin", "chrysobactin", "ornibactin", "precolibactin", "colibactin", "ferripyoverdine"),
                         "pili-related": ("pilin", "pilus", "prepilin", "fimbrial", "fimbriae", "fimbrillin", "fimbrin"),
                         "flagella-related": ("motor", "flagellin", "flagellar", "flagella", "flagellum"),
                         "regulation": ("regulator", "regulatory", "regulation", "receptor", "effector"),
                         "transport": ("transport", "transporter", "permease", "porin", "export"),
                         "toxin": ("toxin", "enterotoxin", "cytotoxic", "necrotizing", "endotoxin", "leukotoxin", "exotoxin", "cytolethal", "o-antigen", "lipooligosaccharide", "lipopolysaccharide"),
                         "chaperones": ("chaperone", "chaperonin"),
                         "cytolysins": ("hemolysin", "cytolysin"),
                         "cell wall related": ("cell wall", ),
                         "cell membrane related": ("membrane", ),
                         "capsule-related": ("capsular", "capsule"),
                         "lipoglycans": ("o-antigen", "lipooligosaccharide", "lipopolysaccharide"),
                         "secretion": ("secretion", "secreted", "secretory"),
                         "efflux": ("efflux", )}

# Counting most common genera
annotation_counter_col_name = "host_genera"
host_genera_counters = CounterWrapper.count_words_in_series(significant_pvals_df[annotation_counter_col_name])
CounterWrapper.dump_counter(host_genera_counters, file="/data1/bio/projects/ndanilova/colitis_crohn/VFDB/pvals/{a}/{b}_significant_{c}_counter.tsv".format(a=pivot_value_col_name, b=group_digest, c=annotation_counter_col_name))

# Dynamic associations dict
host_genera_dict = {i[0].capitalize(): (i[0],) for i in host_genera_counters.most_common(len(gene_description_dict))}

group_data_df = ""
