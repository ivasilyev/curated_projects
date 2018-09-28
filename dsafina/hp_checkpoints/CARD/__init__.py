#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from meta.scripts.Utilities import Utilities
from dsafina.hp_checkpoints.ProjectDescriber import ProjectDescriber
from meta.scripts.card.ReferenceDescriber import ReferenceDescriber
from meta.scripts.LaunchGuideLiner import LaunchGuideLiner

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3
"""

# Prepare new raw reads for filtering to HG19 DB
projectDescriber = ProjectDescriber()
projectDescriber.sampledata = "/data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints_no_hg19.sampledata"
raw_reads_dirs = ["/data1/bio/170922/", "/data1/bio/180419_NB501097_0017_AHL55TBGX5/HP/"]

Utilities.create_sampledata(raw_reads_dirs, projectDescriber.sampledata)
subprocess.getoutput("sed -i '/fastq$/d' {}".format(projectDescriber.sampledata))  # The folders only contains .fastq.gz

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

# Checkout (from WORKER node)
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG python3 \
/home/docker/scripts/verify_coverages.py -s /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints_no_hg19.sampledata \
-r /data/reference/homo_sapiens/hg/hg19/hg19.refdata \
-m hg19 -d -o /data2/bio/Metagenomes/HG19/
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
