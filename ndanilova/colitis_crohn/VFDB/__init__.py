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
from meta.scripts.GroupDataPreparer import GroupDataPreparer
from meta.scripts.RemoteScript import RemoteScript
import subprocess
import os
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
# Charts directory: '/data1/bio/projects/ndanilova/colitis_crohn/VFDB/'

# Look for Redis pod & service:
kubectl get pods

# (Optional) If required, flush all Redis queues:
export IMG=ivasilyev/bwt_filtering_pipeline_master:latest && docker pull $IMG && docker run --rm -it $IMG redis-cli -h redis flushall

# (Optional) Or remove Redis server by itself:
kubectl delete svc,pod redis

# (Optional) Deploy if not present:
kubectl create -f  https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/meta/charts/redis.yaml

# Deploy master chart
kubectl create -f /data1/bio/projects/ndanilova/colitis_crohn/VFDB/master.yaml

# Deploy worker chart only when master pod is complete
kubectl create -f /data1/bio/projects/ndanilova/colitis_crohn/VFDB/worker.yaml

# View active nodes
kubectl describe pod ndanilova-vfdb-job- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/Toxins/VFDB/Statistics/*coverage.txt | wc -l) OF $(cat /data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.groupdata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod ndanilova-vfdb-queue
kubectl delete job ndanilova-vfdb-job

# Checkout (from WORKER node)
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && docker pull $IMG && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG python3 /home/docker/scripts/verify_coverages.py -s /data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.groupdata -r /data/reference/VFDB/vfdb_v2018.06.29/index/vfdb_v2018.06.29.json -m no_hg19 -d -o /data2/bio/Metagenomes/Toxins/VFDB/
"""


groupdata = GroupDataPreparer(raw_groupdata_file=ProjectDescriber.groupdata,
                              processed_prefix="/data2/bio/Metagenomes/Toxins/VFDB/Statistics/",
                              processed_suffix="_no_hg19_vfdb_v2018.06.29_coverage.tsv")

groupdata.export_groupdata("{a}{b}.groupdata".format(a=outputDir,
                                                         b=groupdata.groups_string))

merger = RemoteScript(url="https://raw.githubusercontent.com/ivasilyev/statistical_tools/master/groupdata2statistics.py")
merger.download()

subprocess.getoutput("python3 {a} -g /data1/bio/projects/ndanilova/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.groupdata -i reference_id -v id_mapped_relative_abundance -o {b}pvals".format(a=merger.file, b=outputDir))

annotator = PivotTableAnnotator(pivot_file="/data1/bio/projects/ndanilova/colitis_crohn/VFDB/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr_total_dataframe.tsv",
                                annotation_file="/data/reference/VFDB/vfdb_v2018.06.29/index/vfdb_v2018.06.29_annotation.tsv")

annotator.export_annotated_pivot("/data1/bio/projects/ndanilova/colitis_crohn/VFDB/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr_total_dataframe_annotated.tsv")
annotator.export_group_reports("/data1/bio/projects/ndanilova/colitis_crohn/VFDB/group_reports")
