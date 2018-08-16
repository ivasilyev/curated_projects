#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3
"""

from tgrigoreva.ecoli_plates_suspension.ProjectDescriber import ProjectDescriber
from meta.scripts.vfdb.ReferenceDescriber import ReferenceDescriber
from meta.scripts.LaunchGuideLiner import LaunchGuideLiner
import os

outputDir = "/data1/bio/projects/tgrigoreva/ecoli_plates_suspension/VFDB/"
os.makedirs(outputDir, exist_ok=True)

guideLiner = LaunchGuideLiner(charts_dir="/data1/bio/projects/tgrigoreva/ecoli_plates_suspension/VFDB/charts/",
                              deploy_prefix="tgrigoreva-vfdb",
                              nodes_number=7,
                              threads_number="half",
                              sampledata_file=ProjectDescriber.sampledata,
                              refdata_file=ReferenceDescriber.refdata,
                              output_mask="no_hg19",
                              output_dir="/data2/bio/Metagenomes/Toxins/VFDB/")

guideLiner.generate_config()
guideLiner.get_deploy_guide()

"""
# Charts directory: '/data1/bio/projects/tgrigoreva/ecoli_plates_suspension/VFDB/charts/'

# Look for Redis pod & service:
kubectl get pods

# (Optional) If required, flush all Redis queues:
export IMG=ivasilyev/bwt_filtering_pipeline_master:latest && docker pull $IMG && docker run --rm -it $IMG redis-cli -h redis flushall

# (Optional) Or remove Redis server by itself:
kubectl delete svc,pod redis

# (Optional) Deploy if not present:
kubectl create -f  https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/meta/charts/redis.yaml

# Deploy master chart
kubectl create -f /data1/bio/projects/tgrigoreva/ecoli_plates_suspension/VFDB/charts/master.yaml

# Deploy worker chart only when the master pod ('tgrigoreva-vfdb-queue') is complete
kubectl create -f /data1/bio/projects/tgrigoreva/ecoli_plates_suspension/VFDB/charts/worker.yaml

# View active nodes
kubectl describe pod tgrigoreva-vfdb-job- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/Toxins/VFDB/Statistics/*coverage.txt | wc -l) OF $(cat /data1/bio/projects/tgrigoreva/ecoli_plates_suspension/raw.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod tgrigoreva-vfdb-queue
kubectl delete job tgrigoreva-vfdb-job

# Checkout (from WORKER node)
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && docker pull $IMG && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG python3 /home/docker/scripts/verify_coverages.py -s /data1/bio/projects/tgrigoreva/ecoli_plates_suspension/raw.sampledata -r /data/reference/VFDB/vfdb_v2018.06.29/index/vfdb_v2018.06.29.json -m no_hg19 -d -o /data2/bio/Metagenomes/Toxins/VFDB/

"""
