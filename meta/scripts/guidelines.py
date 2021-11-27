#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import yaml
from meta.scripts.ChartGenerator import ChartGenerator
from meta.scripts.Utilities import Utilities


def dump_index_guide(input_nucleotide_fasta: str, output_dir: str):
    if not Utilities.is_file_valid(input_nucleotide_fasta):
        raise ValueError(f"Invalid file: '{input_nucleotide_fasta}'")
    cmd_0 = f"""
    export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
    docker pull "$IMG" && \
    docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it "$IMG" \
        bash -c '
            cd "{output_dir}";
            python3 "$HOME/scripts/cook_the_reference.py" \
                --input "{input_nucleotide_fasta}" \
                --output "{output_dir}";
        '
    """
    cmd = Utilities.join_lines(cmd_0)
    out_file = os.path.join(output_dir, "index.sh")
    Utilities.dump_string(cmd, out_file)
    print(f"For indexing, run outside of Docker: 'bash \"{out_file}\"'")


class Chart(object):
    def __init__(self, file, url):
        self.file = file
        self.url = url

    def finalize(self, config_file):
        g = ChartGenerator(config=config_file, template=self.url)
        g.export_yaml(self.file)
        print("Created chart: '{}'".format(self.file))


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
        self.charts_directory = Utilities.ends_with_slash(charts_dir)
        self.deploy_prefix = re.sub("[^A-Za-z0-9\-]+", "-", deploy_prefix)
        self.config_chart = Chart(file="{}config.yaml".format(self.charts_directory),
                                  # URL is not supported
                                  url="https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/config.yaml")
        self.cfgDict = {"QUEUE_NAME": "{}-queue".format(self.deploy_prefix),
                        "MASTER_CONTAINER_NAME": "{}-master".format(self.deploy_prefix),
                        "JOB_NAME": "{}-job".format(self.deploy_prefix),
                        "WORKER_CONTAINER_NAME": "{}-worker".format(self.deploy_prefix),
                        "ACTIVE_NODES_NUMBER": nodes_number,
                        "THREADS_NUMBER": threads_number,
                        "SAMPLEDATA": sampledata_file,
                        "REFDATA": refdata_file,
                        "OUTPUT_MASK": output_mask,
                        "OUTPUT_DIR": Utilities.ends_with_slash(output_dir)}
        self.master_chart = Chart(file="{}master.yaml".format(self.charts_directory),
                                  url="https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/master.yaml")
        self.worker_chart = Chart(file="{}worker.yaml".format(self.charts_directory),
                                  url="https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/bwt_filtering_pipeline/templates/bwt-fp-only-coverage/worker.yaml")

    def generate_config(self):
        os.makedirs(self.charts_directory, exist_ok=True)
        cfg = self.config_chart.file
        with open(cfg, mode="w", encoding="utf-8") as f:
            yaml.dump(self.cfgDict, f, default_flow_style=False, explicit_start=False)
            print("Created project configuration: '{}'".format(self.cfgDict))
        self.master_chart.finalize(cfg)
        self.worker_chart.finalize(cfg)

    def get_deploy_guide(self):
        print("""
# Charts directory: '{chartsDir}' 

# Look for Redis pod & service:
kubectl get pods

# (Optional) If required, flush all Redis queues:
export IMG=ivasilyev/bwt_filtering_pipeline_master:latest && \\
docker pull $IMG && docker run --rm -it $IMG redis-cli -h redis flushall

# (Optional) Or remove Redis server by itself:
kubectl delete svc,pod redis

# (Optional) Deploy if not present:
kubectl create -f  https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/meta/charts/redis.yaml

# Deploy master chart
kubectl create -f {MASTER}

# Deploy worker chart only when the master pod ('{QUEUE_NAME}') is complete
kubectl create -f {WORKER}

# View active nodes
kubectl describe pod {JOB_NAME}- | grep Node:

# View progress (from WORKER node)
echo && echo PROCESSED $(ls -d {OUTPUT_DIR}Statistics/*coverage.txt | wc -l) OF $(cat {sampleDataFileName} | wc -l)

# Look for some pod (from MASTER node)
kubectl describe pod <NAME>

# Cleanup
kubectl delete pod {QUEUE_NAME}
kubectl delete job {JOB_NAME}

# Checkout (from WORKER node)
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \\
docker pull $IMG && \\
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \\ 
python3 /home/docker/scripts/verify_coverages.py \\
-i {sampleDataFileName} \\
-p {OUTPUT_DIR}Statistics/ \\
-s _{OUTPUT_MASK}_coverage.tsv \\
-g <GENOME> \\
-d
        """.format(MASTER=self.master_chart.file,
                   WORKER=self.worker_chart.file,
                   chartsDir=self.charts_directory,
                   JOB_NAME=self.cfgDict["JOB_NAME"],
                   OUTPUT_DIR=self.cfgDict["OUTPUT_DIR"],
                   QUEUE_NAME=self.cfgDict["QUEUE_NAME"],
                   sampleDataFileName=self.cfgDict["SAMPLEDATA"],
                   REFDATA=self.cfgDict["REFDATA"],
                   OUTPUT_MASK=self.cfgDict["OUTPUT_MASK"]))
