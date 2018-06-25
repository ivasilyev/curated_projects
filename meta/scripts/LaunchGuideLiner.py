#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import yaml
from meta.scripts.RemoteScript import RemoteScript
from meta.scripts.KubernetesJobChartsGenerator import KubernetesJobChartsGenerator


class LaunchGuideLiner:
    def __init__(self,
                 reference_database_name,
                 reference_index_directory,
                 reference_nfasta_file,
                 project_owner,
                 project_name,
                 charts_directory,
                 deploy_prefix,
                 output_directory,
                 output_mask,
                 sampledata_file,
                 threads,
                 nodes_number):
        self.reference_database_name = reference_database_name
        self.reference_index_directory = self.ends_with_slash(reference_index_directory)
        self.reference_nfasta_file = reference_nfasta_file
        self.project_owner = project_owner
        self.project_name = project_name
        self.charts_directory = self.ends_with_slash(charts_directory)
        self.deploy_prefix = deploy_prefix
        self.output_directory = "{a}{b}".format(a=self.ends_with_slash(output_directory),
                                                b=self.ends_with_slash(self.reference_database_name))
        self.output_mask = output_mask
        self.sampledata_file = sampledata_file
        self.threads = threads
        self.nodes_number = nodes_number
        # Latest config file: https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/config.yaml
        self.cfgDict = {"QUEUE_NAME": "{}-queue".format(self.deploy_prefix),
                        "MASTER_CONTAINER_NAME": "{}-master".format(self.deploy_prefix),
                        "JOB_NAME": "{}-job".format(self.deploy_prefix),
                        "ACTIVE_NODES_NUMBER": self.nodes_number,
                        "THREADS_NUMBER": self.threads,
                        "WORKER_CONTAINER_NAME": "{}-worker".format(self.deploy_prefix),
                        "SAMPLEDATA": self.sampledata_file,
                        "REFDATA": "{a}{b}.refdata".format(a=self.reference_index_directory, b=self.reference_database_name),
                        "OUTPUT_MASK": self.output_mask,
                        "OUTPUT_DIR": self.output_directory}
        self.cfg_file = "{}config.yaml".format(self.charts_directory)
    @staticmethod
    def ends_with_slash(string):
        if string.endswith("/"):
            return string
        else:
            return str(string + "/")
    def get_index_guide(self):
        print("""
              # Reference indexing (from worker node)
              rm -rf {a}
              docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && \
              docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker:latest \
              python3 /home/docker/scripts/cook_the_reference.py \
              -i {b} \
              -o {a}
              """.format(a=self.reference_index_directory, b=self.reference_nfasta_file))
    def generate_config(self):
        os.makedirs(self.charts_directory, exist_ok=True)
        with open(file=self.cfg_file, mode="w", encoding="utf-8") as f:
            yaml.dump(self.cfgDict, f, default_flow_style=False, explicit_start=False)
        gen = KubernetesJobChartsGenerator(config=self.cfg_file,
                                           master="https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/master.yaml",
                                           worker="https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/worker.yaml")
        gen.export_yaml(output_dir=self.charts_directory)
    def get_deploy_guide(self):
        print("""
        # Copy the directory '{chartsDir}' into '{ownerName}/{projectName}' directory and push updates for remote deployment
    
        # Look for Redis pod & service:
        kubectl get pods --show-all
    
        # Deploy if not present:
        kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-pod.yaml && \
        kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-service.yaml
    
        # Pipeline launch
        # Deploy the MASTER chart to create queue
        ## Remote:
        kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/{ownerName}/{projectName}/charts/master.yaml
        ## Local:
        kubectl create -f {chartsDir}master.yaml
            
        # Wait until master finish and deploy the WORKER chart to create the pipeline job
        ## Remote:
        kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/{ownerName}/{projectName}/charts/worker.yaml
        ## Local:
        kubectl create -f {chartsDir}worker.yaml
        
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
        -m {OUTPUT_MASK}_{referenceDBName} -d -o {OUTPUT_DIR}    
        """.format(ownerName=self.project_owner,
                   chartsDir=self.charts_directory,
                   projectName=self.project_name,
                   JOB_NAME=self.cfgDict["JOB_NAME"],
                   OUTPUT_DIR=self.cfgDict["OUTPUT_DIR"],
                   QUEUE_NAME=self.cfgDict["QUEUE_NAME"],
                   sampleDataFileName=self.cfgDict["SAMPLEDATA"],
                   REFDATA=self.cfgDict["REFDATA"],
                   OUTPUT_MASK=self.cfgDict["OUTPUT_MASK"],
                   referenceDBName=self.reference_database_name))
