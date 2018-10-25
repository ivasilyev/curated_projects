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

# Prepare new raw reads for filtering to HG19 DB
projectDescriber = ProjectDescriber()
projectDescriber.sampledata = "/data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints_no_hg19.sampledata"
raw_reads_dirs = ["/data1/bio/170922/",
                  "/data1/bio/180419_NB501097_0017_AHL55TBGX5/HP/",
                  "/data1/bio/181017_NB501097_0024_AHL553BGX5/Conversion/"]

# Create sampledata from Illumina raw reads
raw_reads_dict = {}
for raw_reads_dir in raw_reads_dirs:
    raw_reads_dir = Utilities.ends_with_slash(raw_reads_dir)
    files_list = os.listdir(raw_reads_dir)
    for file_name in files_list:
        if any([file_name.endswith(i) for i in ["csfasta", "fasta", "fa", "fastq", "fq", "gz"]]):
            sample_name = file_name.split("_")[0].strip()
            if "HP" not in sample_name or os.path.isfile("/data2/bio/Metagenomes/HG19/Unmapped_reads/{}_no_hg19.1.gz".format(sample_name)):
                continue
            sample_extension = file_name.split(".")[-1]
            sample_files = sorted(sorted([raw_reads_dir + i for i in files_list if len(re.findall("^{}".format(sample_name), i)) > 0 and i.endswith(sample_extension)], key=len, reverse=True)[:2])
            sample_name = re.sub("[^A-Za-z0-9]+", "_", sample_name)
            sample_name = re.sub("_+", "_", sample_name)
            raw_reads_dict[sample_name] = sample_files

raw_reads_dict = dict(sorted(raw_reads_dict.items()))
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

# (Optional) If required, flush all Redis queues:
export IMG=ivasilyev/bwt_filtering_pipeline_master:latest && \
docker pull $IMG && docker run --rm -it $IMG redis-cli -h redis flushall

# (Optional) Or remove Redis server by itself:
kubectl delete svc,pod redis

# (Optional) Deploy if not present:
kubectl create -f  https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/meta/charts/redis.yaml

# Deploy master chart
kubectl create -f /data1/bio/projects/dsafina/hp_checkpoints/hg19/charts/master.yaml

# Deploy worker chart only when the master pod ('dsafina-hg19-queue') is complete
kubectl create -f /data1/bio/projects/dsafina/hp_checkpoints/hg19/charts/worker.yaml

# View active nodes
kubectl describe pod dsafina-hg19-job- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/HG19/Statistics/*coverage.txt | wc -l) OF $(cat /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints_no_hg19.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod dsafina-hg19-queue
kubectl delete job dsafina-hg19-job
"""

subprocess.getoutput("sed -i 's|/home/docker/scripts/master.py, |/home/docker/scripts/master.py, -n, |g' /data1/bio/projects/dsafina/hp_checkpoints/hg19/charts/master.yaml")

# Deploy from MASTER

# Prepare sample data for CARD DB
projectDescriber = ProjectDescriber()
# projectDescriber.sampledata = "/data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints_new.sampledata"
# raw_reads_dirs = ["/data2/bio/Metagenomes/HG19/Non-mapped_reads/", "/data2/bio/Metagenomes/HG19/Unmapped_reads/"]
# Utilities.create_sampledata(raw_reads_dirs, projectDescriber.sampledata)

referenceDescriber = ReferenceDescriber()

launchGuideLiner = LaunchGuideLiner(charts_dir="{}{}/charts/".format(Utilities.ends_with_slash(projectDescriber.directory), referenceDescriber.alias),
                                    deploy_prefix="{a}-{b}-{c}".format(a=projectDescriber.owner, b=projectDescriber.name, c=referenceDescriber.alias),
                                    nodes_number=7,
                                    threads_number="half",
                                    sampledata_file=projectDescriber.sampledata,
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

# Deploy worker chart only when the master pod ('dsafina-card-v2-0-3-queue') is complete
kubectl create -f /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/charts/worker.yaml

# View active nodes
kubectl describe pod dsafina-card-v2-0-3-job- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/CARD/Statistics/*coverage.txt | wc -l) OF $(cat /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod dsafina-card-v2-0-3-queue
kubectl delete job dsafina-card-v2-0-3-job

# Checkout (from WORKER node)
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/verify_coverages.py \
-i /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints.sampledata \
-p /data2/bio/Metagenomes/CARD/Statistics/ \
-s _card_v2.0.3_coverage.tsv \
-g /data/reference/CARD/card_v2.0.3/index/card_v2.0.3_samtools.genome \
-d
"""

launchGuideLiner = LaunchGuideLiner(charts_dir="{}{}/charts/".format(Utilities.ends_with_slash(projectDescriber.directory), referenceDescriber.alias),
                                    deploy_prefix="{a}-{b}-{c}".format(a=projectDescriber.owner, b=projectDescriber.name, c=referenceDescriber.alias),
                                    nodes_number=7,
                                    threads_number="half",
                                    sampledata_file="/data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-09-27-11-51-43.sampledata",
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
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/CARD/Statistics/*coverage.txt | wc -l) OF $(cat /data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-09-27-11-51-43.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod dsafina-hp-checkpoints-card-v2-0-3-queue
kubectl delete job dsafina-hp-checkpoints-card-v2-0-3-job
"""
launchGuideLiner = LaunchGuideLiner(charts_dir="{}{}/charts/".format(Utilities.ends_with_slash(projectDescriber.directory), referenceDescriber.alias),
                                    deploy_prefix="{a}-{b}-{c}".format(a=projectDescriber.owner, b=projectDescriber.name, c=referenceDescriber.alias),
                                    nodes_number=7,
                                    threads_number="half",
                                    sampledata_file="/data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-09-27-21-50-35.sampledata",
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
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/CARD/Statistics/*coverage.txt | wc -l) OF $(cat /data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-09-27-21-50-35.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod dsafina-hp-checkpoints-card-v2-0-3-queue
kubectl delete job dsafina-hp-checkpoints-card-v2-0-3-job

"""

launchGuideLiner = LaunchGuideLiner(charts_dir="{}{}/charts/".format(Utilities.ends_with_slash(projectDescriber.directory), referenceDescriber.alias),
                                    deploy_prefix="{a}-{b}-{c}".format(a=projectDescriber.owner, b=projectDescriber.name, c=referenceDescriber.alias),
                                    nodes_number=7,
                                    threads_number="half",
                                    sampledata_file="/data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-09-28-09-21-09.sampledata",
                                    refdata_file=referenceDescriber.refdata,
                                    output_mask=referenceDescriber.alias,
                                    output_dir="/data2/bio/Metagenomes/{}/".format(referenceDescriber.name))
launchGuideLiner.generate_config()
launchGuideLiner.get_deploy_guide()

# Rerun with charts from '2018-09-28-09-21-09.sampledata'
# Again, regenerate charts & manual launch for '2018-09-28-11-12-45.sampledata':
"""
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/nBee.py \
-i /data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-09-28-11-12-45.sampledata \
-r /data/reference/CARD/card_v2.0.3/index/card_v2.0.3_refdata.json \
-m card_v2.0.3 -t half -o /data2/bio/Metagenomes/CARD/
"""

"""
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/nBee.py \
-i /data2/bio/Metagenomes/CARD/Statistics/sampledata/2018-09-29-02-42-19.sampledata \
-r /data/reference/CARD/card_v2.0.3/index/card_v2.0.3_refdata.json \
-m card_v2.0.3 -t half -o /data2/bio/Metagenomes/CARD/
"""

"""
# Combine data for RPM

export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git clone https://github.com/ivasilyev/statistical_tools.git
cd statistical_tools

python3 groupdata2statistics.py -g /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints_new.groupdata \
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

python3 groupdata2statistics.py -g /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints_new.groupdata \
-p /data2/bio/Metagenomes/CARD/Statistics/ \
-s _card_v2.0.3_coverage.tsv \
-i reference_id \
-v id_mapped_reads_per_kbp_per_million_sample_mapped_reads \
-o /data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/pvals/RPKM/
"""

import pandas as pd

index_col_name = "reference_id"

annotation_df = pd.read_table("/data/reference/CARD/card_v2.0.3/index/card_v2.0.3_annotation.tsv").set_index(index_col_name)
for total_df_file in ["/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/pvals/RPM/1_2_3_C_srr_total_dataframe.tsv",
                      "/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/pvals/RPKM/1_2_3_C_srr_total_dataframe.tsv"]:
    total_df = pd.read_table(total_df_file).set_index(index_col_name)
    annotated_total_df = pd.concat([annotation_df, total_df], axis=1)
    annotated_total_df.index.name = index_col_name
    annotated_total_df.to_csv(total_df_file.replace("_total_dataframe.tsv", "_total_dataframe_annotated.tsv"), sep='\t', header=True, index=True)

# TODO Iterate it
annotated_total_df = pd.read_table("/data1/bio/projects/dsafina/hp_checkpoints/card_v2.0.3/pvals/RPM/1_2_3_C_srr_total_dataframe_annotated.tsv").set_index(index_col_name)


class ReversedGroupComparator:
    def __init__(self, name):
        self.name = str(name)
        self.columns = []
        self.rows = []
        self.comparisons = []


group_names = ["C", "1", "2", "3"]
groupdata_comparators_dict = {}
for first_group_name in group_names:
    reversedGroupComparator = ReversedGroupComparator(first_group_name)
    for annotated_total_df_column_name in list(annotated_total_df):
        annotated_total_df_column_name_digest_list = annotated_total_df_column_name.split("_")
        if annotated_total_df_column_name_digest_list[0].strip() == first_group_name and len(annotated_total_df_column_name_digest_list) > 0 and annotated_total_df_column_name_digest_list[1].strip().startswith("/"):
            reversedGroupComparator.columns.append(annotated_total_df_column_name)
    reversedGroupComparator.rows = annotated_total_df.loc[annotated_total_df["{}_non-zero_values".format(first_group_name)].astype(int) > 0, :].index.tolist()
    for second_group_name in [i for i in group_names if i != first_group_name]:
        comparison_name = "{}_vs_{}_is_rejected_by_fdr_bh_for_wilcoxon".format(first_group_name, second_group_name)
        if comparison_name in list(annotated_total_df):
            reversedGroupComparator.comparisons.append(comparison_name)
    groupdata_comparators_dict[first_group_name] = reversedGroupComparator

drug_classes = ["aminoglycoside", "fluoroquinolone", "glycopeptide antibiotic", "lincosamide", "macrolide", "nucleoside antibiotic", "penam", "peptide antibiotic", "phenicol", "sulfonamide", "tetracycline", "triclosan"]
resistance_mechanism = ["efflux", "inactivation", "reduced permeability", "target alteration", "target protection", "target replacement"]

