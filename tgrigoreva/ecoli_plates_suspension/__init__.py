#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${DOCKER_IMAGE_NAME} python3

"""

import os
import subprocess
import yaml
import pandas as pd
import re


subProjectName = "ecoli_plates_suspension"
projectName = "tgrigoreva/{}".format(subProjectName)
outputDir = "/data1/bio/projects/{}/".format(projectName)
os.makedirs(outputDir, exist_ok=True)


def filename_only(string):
    return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))


def list_to_file(header, list_to_write, file_to_write):
    header += "\n".join(str(i) for i in list_to_write if i is not None) + "\n"
    file = open(file_to_write, 'w')
    file.write(header)
    file.close()


# Create sampledata with files of non-zero length only
sampleDataList = ["{}\t{}".format(filename_only(i), i) for i in subprocess.getoutput("ls -d /data1/bio/ecoli_komfi/fasta/*.fasta").split("\n") if len(i) > 0 and int(subprocess.getoutput("wc -l < " + i).split('\n')[0]) > 0]
sampleDataFileName = outputDir + "{}.sampledata".format(subProjectName)
list_to_file("", sampleDataList, sampleDataFileName)


def external_route(*args):
    process = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    print(output.decode("utf-8"))


# Create config chart
chartsSubDir = "charts/"
chartsDir = outputDir + chartsSubDir
os.makedirs(chartsDir, exist_ok=True)
cfgDict = {"QUEUE_NAME": "tgrigoreva-bwt-eps-queue",
           "MASTER_CONTAINER_NAME": "tgrigoreva-bwt-eps-master",
           "JOB_NAME": "tgrigoreva-bwt-eps-job",
           "ACTIVE_NODES_NUMBER": 9,
           "WORKER_CONTAINER_NAME": "tgrigoreva-bwt-eps-worker",
           "SAMPLEDATA": sampleDataFileName,
           "REFDATA": "/data/reference/custom/25_ecoli_genes/index/25_ecoli_genes.refdata",
           "OUTPUT_MASK": "no_hg19",
           "OUTPUT_DIR": "/data2/bio/Metagenomes/custom/25_ecoli_genes"}
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
-m {OUTPUT_MASK}_{subProjectName} -d -o {OUTPUT_DIR}

""".format(chartsDir=chartsDir,
           chartsSubDir=chartsSubDir,
           projectName=projectName,
           JOB_NAME=cfgDict["JOB_NAME"],
           OUTPUT_DIR=cfgDict["OUTPUT_DIR"],
           QUEUE_NAME=cfgDict["QUEUE_NAME"],
           sampleDataFileName=sampleDataFileName,
           REFDATA=cfgDict["REFDATA"],
           OUTPUT_MASK=cfgDict["OUTPUT_MASK"],
           subProjectName=subProjectName))

"""
# Copy the directory '/data1/bio/projects/tgrigoreva/ecoli_plates_suspension/charts/' into 'tgrigoreva/ecoli_plates_suspension' directory and push updates

# Look for Redis pod & service:
kubectl get pods --show-all

# Deploy if not present:
kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-pod.yaml && kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-service.yaml

# Pipeline launch
# Deploy the MASTER chart to create queue
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/tgrigoreva/ecoli_plates_suspension/charts/master.yaml

# Wait until master finish and deploy the WORKER chart to create the pipeline job
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/tgrigoreva/ecoli_plates_suspension/charts/worker.yaml

# View active nodes
kubectl describe pod tgrigoreva-bwt-eps-job- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/custom/25_ecoli_genes/Statistics/*coverage.txt | wc -l) OF $(cat /data1/bio/projects/tgrigoreva/ecoli_plates_suspension/ecoli_plates_suspension.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod tgrigoreva-bwt-eps-queue
kubectl delete job tgrigoreva-bwt-eps-job

# Checkout (from WORKER node)
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker python3 /home/docker/scripts/verify_coverages.py -s /data1/bio/projects/tgrigoreva/ecoli_plates_suspension/ecoli_plates_suspension.sampledata -r /data/reference/custom/25_ecoli_genes/index/25_ecoli_genes.refdata -m no_hg19_ecoli_plates_suspension -d -o /data2/bio/Metagenomes/custom/25_ecoli_genes

"""

sampleDataDF = pd.read_table(sampleDataFileName, sep='\t', header="infer", names=["sample_name", "sample_path"], engine='python')
sampleDataDF["coverage_path"] = sampleDataDF.loc[:, "sample_name"].apply(lambda x: "/data2/bio/Metagenomes/custom/25_ecoli_genes/Statistics/{}_no_hg19_25_ecoli_genes_coverage.txt".format(x))
sampleDataDF = sampleDataDF.loc[sampleDataDF["coverage_path"].map(lambda x: os.path.isfile(x))].set_index("sample_name")

sampleDataDF.index.tolist()
referenceAnnotationDF = pd.read_table("/data/reference/custom/25_ecoli_genes/index/25_ecoli_genes_annotation.txt", sep='\t', header=0, engine='python').set_index("reference_id")
os.makedirs("{}coverages/".format(outputDir), exist_ok=True)


def annotate_df(df, output_file_name):
    df = pd.concat([referenceAnnotationDF, df], axis=1)
    df.to_csv(output_file_name, sep='\t', index=True, header=True)


for ix in sampleDataDF.index.tolist():
    annotate_df(pd.read_table(sampleDataDF.loc[ix, "coverage_path"], sep='\t', header=0, engine='python').set_index("reference_id"), "{}coverages/{}_merged.tsv".format(outputDir, ix))

# Create reads sampledata
sampleDataList = ["\t".join([re.sub("_R1_.*", "_reads", filename_only(i)), i, i.replace("_R1_", "_R2_")]) for i in subprocess.getoutput("ls -d /data1/bio/ecoli_komfi/raw_reads/*.fastq.gz").split("\n") if len(i) > 0 and "_R1_" in i and os.path.isfile(i.replace("_R1_", "_R2_"))]
sampleDataFileName = outputDir + "{}_reads.sampledata".format(subProjectName)
list_to_file("", sampleDataList, sampleDataFileName)

# Create next alignment config chart
chartsSubDir = "reads_charts/"
chartsDir = outputDir + chartsSubDir
os.makedirs(chartsDir, exist_ok=True)
cfgDict = {"QUEUE_NAME": "tgrigoreva-bwt-eps2-queue",
           "MASTER_CONTAINER_NAME": "tgrigoreva-bwt-eps2-master",
           "JOB_NAME": "tgrigoreva-bwt-eps2-job",
           "ACTIVE_NODES_NUMBER": 9,
           "WORKER_CONTAINER_NAME": "tgrigoreva-bwt-ep2s-worker",
           "SAMPLEDATA": sampleDataFileName,
           "REFDATA": "/data/reference/custom/25_ecoli_genes/index/25_ecoli_genes.refdata",
           "OUTPUT_MASK": "no_hg19",
           "OUTPUT_DIR": "/data2/bio/Metagenomes/custom/25_ecoli_genes"}
cfgFileName = chartsDir + "config.yaml"
with open(cfgFileName, 'w') as cfgFile:
    yaml.dump(cfgDict, cfgFile, default_flow_style=False, explicit_start=False)

genFileName = chartsDir + "generator.py"
external_route("curl", "-fsSL", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/generator.py", "-o", genFileName)

# Create charts from templates
external_route("python3", genFileName, "-c", cfgFileName, "-m", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/master.yaml", "-w", "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/worker.yaml", "-o", chartsDir)
os.remove(genFileName)

"""
# Copy the directory '/data1/bio/projects/tgrigoreva/ecoli_plates_suspension/reads_charts/' into 'tgrigoreva/ecoli_plates_suspension' directory and push updates

# Look for Redis pod & service:
kubectl get pods --show-all

# Deploy if not present:
kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-pod.yaml && kubectl create -f https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/test_charts/redis-service.yaml

# Pipeline launch
# Deploy the MASTER chart to create queue
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/tgrigoreva/ecoli_plates_suspension/reads_charts/master.yaml

# Wait until master finish and deploy the WORKER chart to create the pipeline job
kubectl create -f https://raw.githubusercontent.com/ivasilyev/curated_projects/master/tgrigoreva/ecoli_plates_suspension/reads_charts/worker.yaml

# View active nodes
kubectl describe pod tgrigoreva-bwt-eps2-job- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d /data2/bio/Metagenomes/custom/25_ecoli_genes/Statistics/*coverage.txt | wc -l) OF $(cat /data1/bio/projects/tgrigoreva/ecoli_plates_suspension/ecoli_plates_suspension_reads.sampledata | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod tgrigoreva-bwt-eps2-queue
kubectl delete job tgrigoreva-bwt-eps2-job

# Checkout (from WORKER node)
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker python3 /home/docker/scripts/verify_coverages.py -s /data1/bio/projects/tgrigoreva/ecoli_plates_suspension/ecoli_plates_suspension_reads.sampledata -r /data/reference/custom/25_ecoli_genes/index/25_ecoli_genes.refdata -m no_hg19_ecoli_plates_suspension -d -o /data2/bio/Metagenomes/custom/25_ecoli_genes

"""

sampleDataDF = pd.read_table(sampleDataFileName, sep='\t', header="infer", names=["sample_name", "sample_path_1", "sample_path_2"], engine='python')
sampleDataDF["coverage_path"] = sampleDataDF.loc[:, "sample_name"].apply(lambda x: "/data2/bio/Metagenomes/custom/25_ecoli_genes/Statistics/{}_no_hg19_25_ecoli_genes_coverage.txt".format(x))
sampleDataDF = sampleDataDF.loc[sampleDataDF["coverage_path"].map(lambda x: os.path.isfile(x))].set_index("sample_name")

for ix in sampleDataDF.index.tolist():
    annotate_df(pd.read_table(sampleDataDF.loc[ix, "coverage_path"], sep='\t', header=0, engine='python').set_index("reference_id"), "{}coverages/{}_annotated.tsv".format(outputDir, ix))
