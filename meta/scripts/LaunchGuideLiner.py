#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import yaml
from meta.scripts.ChartGenerator import ChartGenerator
from meta.scripts.utilities import ends_with_slash


class LaunchGuideLiner:
    def __init__(self,
                 charts_dir,
                 deploy_prefix,
                 nodes_number,
                 threads_number,
                 sampledata_file,
                 refdata_file,
                 output_mask,
                 output_dir):
        self.charts_directory = ends_with_slash(charts_dir)
        self.deploy_prefix = deploy_prefix
        # Latest config sample: https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/config.yaml
        self.cfgDict = {"QUEUE_NAME": "{}-queue".format(deploy_prefix),
                        "MASTER_CONTAINER_NAME": "{}-master".format(deploy_prefix),
                        "JOB_NAME": "{}-job".format(deploy_prefix),
                        "WORKER_CONTAINER_NAME": "{}-worker".format(deploy_prefix),
                        "ACTIVE_NODES_NUMBER": nodes_number,
                        "THREADS_NUMBER": threads_number,
                        "SAMPLEDATA": sampledata_file,
                        "REFDATA": refdata_file,
                        "OUTPUT_MASK": output_mask,
                        "OUTPUT_DIR": ends_with_slash(output_dir)}
        self.cfg_file = "{}config.yaml".format(self.charts_directory)
    @staticmethod
    def get_index_guide(index_directory, raw_nfasta_file):
        print("""
# Reference indexing (from worker node):

rm -rf {a}
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker:latest \
python3 /home/docker/scripts/cook_the_reference.py \
-i {b} \
-o {a}

Wait until REFDATA file creates
              """.format(a=index_directory, b=raw_nfasta_file))
    def generate_config(self):
        os.makedirs(self.charts_directory, exist_ok=True)
        with open(file=self.cfg_file, mode="w", encoding="utf-8") as f:
            yaml.dump(self.cfgDict, f, default_flow_style=False, explicit_start=False)
        gen = ChartGenerator(config=self.cfg_file,
                             template="https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/master_worker.yaml")
        gen.export_yaml("{}master_worker.yaml".format(self.charts_directory))
    def get_deploy_guide(self):
        print("""
# Charts directory: '{chartsDir}' 

# Look for Redis pod & service:
kubectl get pods --show-all

# Deploy if not present:
https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/meta/charts/redis.yaml
# Pipeline launch
# Deploy the combined chart to create queue and start pipeline job
kubectl create -f {chartsDir}master_worker.yaml

# View active nodes
kubectl describe pod {JOB_NAME}- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d {OUTPUT_DIR}/Statistics/*coverage.txt | wc -l) OF $(cat {sampleDataFileName} | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod {QUEUE_NAME}
kubectl delete job {JOB_NAME}

# Checkout (from WORKER node)
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker python3 \
/home/docker/scripts/verify_coverages.py -s {sampleDataFileName} \
-r {REFDATA} \
-m {OUTPUT_MASK} -d -o {OUTPUT_DIR}
        """.format(chartsDir=self.charts_directory,
                   JOB_NAME=self.cfgDict["JOB_NAME"],
                   OUTPUT_DIR=self.cfgDict["OUTPUT_DIR"],
                   QUEUE_NAME=self.cfgDict["QUEUE_NAME"],
                   sampleDataFileName=self.cfgDict["SAMPLEDATA"],
                   REFDATA=self.cfgDict["REFDATA"],
                   OUTPUT_MASK=self.cfgDict["OUTPUT_MASK"]))
