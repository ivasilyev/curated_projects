#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${DOCKER_IMAGE_NAME} python3
"""

import subprocess
import os
import yaml
import git


"""
Download data via chinese proxy from http://www.mgc.ac.cn/VFs/download.htm

http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz
http://www.mgc.ac.cn/VFs/Down/Comparative_tables_from_VFDB.tar.gz
http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz
http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz

"""

ownerName = "dsafina"
projectName = "hp_checkpoints"

# VFDB updates at every Friday
def get_last_friday():
    import datetime
    import calendar
    last_friday = datetime.date.today()
    oneday = datetime.timedelta(days=1)
    while last_friday.weekday() != calendar.FRIDAY:
        last_friday -= oneday
    return "{:%Y.%m.%d}".format(last_friday)


referenceDBName = "vfdb_v{}".format(get_last_friday())
outputDir = "/data1/bio/projects/{a}/{b}/{c}/".format(a=ownerName, b=projectName, c=referenceDBName)
referenceDir = "/data/reference/VFDB/"
referenceSequenceFile = "{referenceDir}{referenceDBName}.fasta".format(referenceDir=referenceDir, referenceDBName=referenceDBName)
subprocess.getoutput("mkdir -p {} {}".format(outputDir, referenceDir))

# Get reference
subprocess.getoutput("curl -s http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz | zcat > {}".format(referenceSequenceFile))

print("""
# Reference indexing
rm -rf {referenceDir}index
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker:latest \
python3 /home/docker/scripts/cook_the_reference.py \
-i {referenceSequenceFile} \
-o {referenceDir}index
""".format(referenceDir=referenceDir, referenceSequenceFile=referenceSequenceFile))

"""
# Reference indexing
rm -rf /data/reference/VFDB/index
docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker:latest python3 /home/docker/scripts/cook_the_reference.py -i /data/reference/VFDB/vfdb_v2018.05.25.fasta -o /data/reference/VFDB/index

"""

chartsDir = "{}charts/".format(outputDir)
deployName = "dsafina-bwt-vfdb"

cfgDict = {"QUEUE_NAME": "{}-queue".format(deployName),
           "MASTER_CONTAINER_NAME": "{}-master".format(deployName),
           "JOB_NAME": "{}-job".format(deployName),
           "ACTIVE_NODES_NUMBER": 9,
           "WORKER_CONTAINER_NAME": "{}-worker".format(deployName),
           "SAMPLEDATA": "/data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints.sampledata",
           "REFDATA": "{referenceDir}index/{referenceDBName}.refdata".format(referenceDir=referenceDir, referenceDBName=referenceDBName),
           "OUTPUT_MASK": "no_hg19",
           "OUTPUT_DIR": "/data2/bio/Metagenomes/Toxins/TADB"}

cfgFileName = "{}config.yaml".format(chartsDir)

subprocess.getoutput("mkdir -p {}".format(chartsDir))
with open(cfgFileName, 'w') as cfgFile:
    yaml.dump(cfgDict, cfgFile, default_flow_style=False, explicit_start=False)

genFileName = chartsDir + "generator.py"
subprocess.getoutput("curl -fsSL https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/generator.py -o {}".format(genFileName))
subprocess.getoutput("python3 {genFileName} -c {cfgFileName} -m https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/master.yaml -w https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/worker.yaml -o {chartsDir}".format(genFileName=genFileName, cfgFileName=cfgFileName, chartsDir=chartsDir))
os.remove(genFileName)
