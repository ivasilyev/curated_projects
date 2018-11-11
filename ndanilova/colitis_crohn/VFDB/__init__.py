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
import subprocess
import os
import re
import pandas as pd
import numpy as np
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

"""
# Combine data for RPM
rm -rf /data1/bio/projects/ndanilova/colitis_crohn/VFDB/pvals /data1/bio/projects/ndanilova/colitis_crohn/VFDB/metadata_digest

export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git clone https://github.com/ivasilyev/statistical_tools.git
cd statistical_tools

python3 groupdata2statistics.py \
-g /data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.groupdata \
-p /data2/bio/Metagenomes/Toxins/VFDB/Statistics/ \
-s _vfdb_v2018.11.09_coverage.tsv \
-i reference_id \
-v id_mapped_reads_per_million_sample_mapped_reads \
-o /data1/bio/projects/ndanilova/colitis_crohn/VFDB/pvals/RPM/
"""


class VfdbHeaderExtractor:
    @staticmethod
    def is_string_valid(s: str):
        s = s.strip()
        if len(s) > 0 and s != "-":
            return s
    @staticmethod
    def process_regex(pattern: str, s: str):
        out = re.findall(pattern, s)
        if len(out) > 0:
            out = VfdbHeaderExtractor.is_string_valid(out[0])
            if out:
                return out
        return np.nan
    @staticmethod
    def get_genbank_id(s: str):
        return VfdbHeaderExtractor.process_regex("^VFG[0-9]+\(([^\(\)]+)\)", s)
    @staticmethod
    def get_gene_name(s: str):
        return VfdbHeaderExtractor.process_regex("^VFG[^ ]+ \(([^\(\)]+)\)", s)
    @staticmethod
    def get_gene_description(s: str):
        return VfdbHeaderExtractor.process_regex("\)([^\[\]\(\)]+)\[", s)


class CounterWrapper:
    @staticmethod
    def prepare_string(s: str):
        return re.sub(" +", " ", re.sub("[^a-z0-9]+", " ", s.lower()))
    @staticmethod
    def count_words_in_series(series: pd.Series):
        lst = [i for i in CounterWrapper.prepare_string(" ".join(series.fillna("").values.tolist())).split(" ") if len(i) > 3]
        return Counter(lst)
    @staticmethod
    def dump_counter(counter: Counter, file: str):
        Utilities.dump_2d_array([("keyword", "occurrences")] + counter.most_common(), file=file)


index_col_name = "reference_id"
annotation_df = pd.read_table("/data/reference/VFDB/vfdb_v2018.11.09/index/vfdb_v2018.11.09_annotation.tsv").set_index(index_col_name).sort_index()
annotation_df["host_strain"] = annotation_df["former_id"].apply(lambda x: re.findall("\[([^\[\]]+)\]$", x.strip())[0].strip())
annotation_df["host_species"] = annotation_df["host_strain"].apply(lambda x: re.findall("([A-Z]{1}[a-z]+ [a-z]+)", x)[0].strip())
annotation_df["host_genera"] = annotation_df["host_strain"].apply(lambda x: re.findall("([A-Z]{1}[a-z]+)", x)[0].strip())
annotation_df["vfdb_id"] = annotation_df["former_id"].apply(lambda x: re.findall("(^VFG[0-9]+)", x)[0].strip())
annotation_df["genbank_id"] = annotation_df["former_id"].apply(VfdbHeaderExtractor.get_genbank_id)
annotation_df["gene_name"] = annotation_df["former_id"].apply(VfdbHeaderExtractor.get_gene_name)
annotation_df["gene_description"] = annotation_df["former_id"].apply(VfdbHeaderExtractor.get_gene_description)
annotation_df = annotation_df.loc[:, [i for i in list(annotation_df) if i != "id_bp"] + ["id_bp"]]

# TODO Place the annotator into the corresponding ReferenceDescriber

pivot_value_col_name = "RPM"
group_digest = "colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr"
total_df_file = "/data1/bio/projects/ndanilova/colitis_crohn/VFDB/pvals/{a}/{b}_total_dataframe.tsv".format(a=pivot_value_col_name, b=group_digest)
total_df = pd.read_table(total_df_file).set_index(index_col_name).sort_index()
annotated_total_df = pd.concat([annotation_df, total_df], axis=1, sort=True)
annotated_total_df.index.name = index_col_name
annotated_total_df.to_csv(total_df_file.replace("_total_dataframe.tsv", "_total_dataframe_annotated.tsv"), sep='\t',
                          header=True, index=True)
significant_pvals_df = annotated_total_df.loc[annotated_total_df["null_hypothesis_rejections_counter"].astype(int) > 0]

annotation_counter_col_name = "gene_description"
gene_description_counters = CounterWrapper.count_words_in_series(significant_pvals_df[annotation_counter_col_name])
CounterWrapper.dump_counter(gene_description_counters, file="/data1/bio/projects/ndanilova/colitis_crohn/VFDB/pvals/{a}/{b}_significant_{c}_counter.tsv".format(a=pivot_value_col_name, b=group_digest, c=annotation_counter_col_name))

re.sub("[^a-z0-9]+", " ", "".lower())

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

annotation_counter_col_name = "host_genera"
host_genera_counters = CounterWrapper.count_words_in_series(significant_pvals_df[annotation_counter_col_name])
CounterWrapper.dump_counter(host_genera_counters, file="/data1/bio/projects/ndanilova/colitis_crohn/VFDB/pvals/{a}/{b}_significant_{c}_counter.tsv".format(a=pivot_value_col_name, b=group_digest, c=annotation_counter_col_name))

host_genera_dict = {i[0].capitalize(): (i[0],) for i in host_genera_counters.most_common(len(gene_description_dict))}


