#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host --env="DISPLAY" -it ${DOCKER_IMAGE_NAME} python3

"""

import re
import multiprocessing
import os
import requests
import bs4
import subprocess
import pandas as pd
import yaml


subProjectName = "SCFAs_from_IGC_by_KEGG"
projectName = "ndanilova/colitis_crohn/{}".format(subProjectName)
outputDir = "/data1/bio/projects/{}/".format(projectName)
os.makedirs(outputDir, exist_ok=True)

# Look "SCFAs_from_KEGG" project for details
groupDataFileName = "/data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.groupdata"
sampleDataFileName = "/data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.sampledata"


"""
# Reference indexing
ln -s /data/reference/IGC/760MetaHit_139HMP_368PKU_511Bac.fa.90_95 /data/reference/IGC/igc_v2014.03.fasta
rm -rf /data/reference/IGC/index
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker:latest \
python3 /home/docker/scripts/cook_the_reference.py \
-i /data/reference/IGC/igc_v2014.03.fasta \
-o /data/reference/IGC/index

"""


def external_route(*args):
    process = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    print(output.decode("utf-8"))


def filename_only(string):
    return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))


# Create config chart
chartsSubDir = "charts/"
chartsDir = outputDir + chartsSubDir
os.makedirs(chartsDir, exist_ok=True)
subProjectAlias = "ndanilova-bwt-igc"
cfgDict = {"QUEUE_NAME": "{}-queue".format(subProjectAlias),
           "MASTER_CONTAINER_NAME": "{}-master".format(subProjectAlias),
           "JOB_NAME": "{}-job".format(subProjectAlias),
           "ACTIVE_NODES_NUMBER": 7,
           "WORKER_CONTAINER_NAME": "{}-worker".format(subProjectAlias),
           "SAMPLEDATA": sampleDataFileName,
           "REFDATA": "/data/reference/IGC/index/igc_v2014.03.refdata",
           "OUTPUT_MASK": "no_hg19",
           "OUTPUT_DIR": "/data2/bio/Metagenomes/IGC"}
referenceFASTAFileName = filename_only(subprocess.getoutput("cat {}".format(cfgDict["REFDATA"])).split("\n")[0].split("\t")[0])
cfgFileName = chartsDir + "config.yaml"
with open(cfgFileName, 'w') as cfgFile:
    yaml.dump(cfgDict, cfgFile, default_flow_style=False, explicit_start=False)

genFileName = chartsDir + "generator.py"
external_route("curl", "-fsSL", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/generator.py", "-o", genFileName)

# Create charts from templates
external_route("python3", genFileName, "-c", cfgFileName, "-m", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/master.yaml", "-w", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/worker.yaml", "-o", chartsDir)
os.remove(genFileName)

print("""
# Copy the directory '{chartsDir}' into '{projectName}' directory and push updates

# Look for Redis pod & service:
kubectl get pods --show-all

# Deploy if not present:
kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-pod.yaml && \
kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-service.yaml

# Pipeline launch
# Deploy the MASTER chart to create queue
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/{projectName}/{chartsSubDir}master.yaml

# Wait until master finish and deploy the WORKER chart to create the pipeline job
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/{projectName}/{chartsSubDir}worker.yaml

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
-m {OUTPUT_MASK}_{referenceFASTAFileName} -d -o {OUTPUT_DIR}

""".format(chartsDir=chartsDir,
           chartsSubDir=chartsSubDir,
           projectName=projectName,
           JOB_NAME=cfgDict["JOB_NAME"],
           OUTPUT_DIR=cfgDict["OUTPUT_DIR"],
           QUEUE_NAME=cfgDict["QUEUE_NAME"],
           sampleDataFileName=sampleDataFileName,
           REFDATA=cfgDict["REFDATA"],
           OUTPUT_MASK=cfgDict["OUTPUT_MASK"],
           referenceFASTAFileName=referenceFASTAFileName))


"""
# Copy the directory '/data1/bio/projects/ndanilova/colitis_crohn/SCFAs_from_IGC_by_KEGG/charts/' into 'ndanilova/colitis_crohn/SCFAs_from_IGC_by_KEGG' directory and push updates

# Look for Redis pod & service:
kubectl get pods --show-all

# Deploy if not present:
kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-pod.yaml && kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-service.yaml

# Pipeline launch
# Deploy the MASTER chart to create queue
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ndanilova/colitis_crohn/SCFAs_from_IGC_by_KEGG/charts/master.yaml

# Wait until master finish and deploy the WORKER chart to create the pipeline job
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ndanilova/colitis_crohn/SCFAs_from_IGC_by_KEGG/charts/worker.yaml

# View active nodes
kubectl describe pod ndanilova-bwt-igc-job- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/IGC/Statistics/*coverage.txt | wc -l) OF $(cat /data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod ndanilova-bwt-igc-queue
kubectl delete job ndanilova-bwt-igc-job

# Checkout (from WORKER node)
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker python3 /home/docker/scripts/verify_coverages.py -s /data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.sampledata -r /data/reference/IGC/index/igc_v2014.03.refdata -m no_hg19_igc_v2014.03 -d -o /data2/bio/Metagenomes/IGC

"""
