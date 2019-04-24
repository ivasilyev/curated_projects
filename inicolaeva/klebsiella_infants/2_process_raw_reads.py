#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# From another console:
export IMG=ivasilyev/spades_cutadapt:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

# Inside of the created container:
cd ${HOME} && \
git clone https://github.com/ivasilyev/curated_projects.git && \
cd ./curated_projects && \
pip3 install pandas && \
python3
"""

import os
import subprocess
import pandas as pd
from meta.scripts.Utilities import Utilities
from inicolaeva.klebsiella_infants.ProjectDescriber import ProjectDescriber


def run_cutadapt(input_list: list):
    # Note the order, it hardly depends from the order of the upstream dataframe columns
    sample_name, sample_file_1, sample_file_2 = input_list
    _ADAPTER = "AGATCGGAAGAG"
    out_file_1, out_file_2, log_file = [os.path.join(cutadaptDir, "{}_cutadapt.{}".format(sample_name, i)) for i in
                                        ("1.fq.gz", "2.fq.gz", "log")]
    cmd = "cutadapt -a {ad} -A {ad} -m 50 -o {o1} -p {o2} {i1} {i2}".format(ad=_ADAPTER, i1=sample_file_1,
                                                                            i2=sample_file_2, o1=out_file_1,
                                                                            o2=out_file_2)
    try:
        for _f in [out_file_1, out_file_2, log_file]:
            if os.path.exists(_f):
                os.remove(_f)
        log = subprocess.getoutput(cmd)
    except PermissionError:
        raise ValueError("Permission denied, please run `sudo chmod -R 777 {}`".format(os.path.dirname(sample_file_1)))
    Utilities.dump_string(log, file=log_file)
    return {"sample_name": sample_name, "trimmed_file_1": out_file_1, "trimmed_file_2": out_file_2}


def run_spades(input_list: list):
    # Same about the order
    sample_name, sample_file_1, sample_file_2 = input_list
    out_dir = os.path.join(spadesDir, sample_name)
    subprocess.getoutput("rm -rf {}".format(out_dir))
    os.makedirs(out_dir)
    cmd = "spades.py --careful -o {out} -1 {i1} -2 {i2}".format(out=out_dir, i1=sample_file_1, i2=sample_file_2)
    log = subprocess.getoutput(cmd)
    log_file = os.path.join(out_dir, "{}_spades.log".format(sample_name))
    Utilities.dump_string(log, file=log_file)
    return {"sample_name": sample_name, "assembly": os.path.join(out_dir, "contigs.fasta")}


projectDescriber = ProjectDescriber()
rawSampledataDF = Utilities.load_tsv("/data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/raw_reads.sampledata")
# Prepare path
rawReadsDir = os.path.join(projectDescriber.RAW_DATA_DIR, "reads")
cutadaptDir = os.path.join(rawReadsDir, "cutadapt")
os.makedirs(cutadaptDir, exist_ok=True)
# Trim reads
cutadaptResults = Utilities.multi_core_queue(run_cutadapt, queue=rawSampledataDF.values.tolist())
cutadaptResultsDF = pd.DataFrame.from_dict(cutadaptResults).sort_values("sample_name")
Utilities.dump_tsv(cutadaptResultsDF, table_file=projectDescriber.SAMPLE_DATA_FILE)
# Assemble reads
spadesDir = os.path.join(rawReadsDir, "spades")
spadesResults = Utilities.single_core_queue(run_spades, cutadaptResultsDF.values.tolist())
spadesResultsSampleData = os.path.join(os.path.dirname(projectDescriber.SAMPLE_DATA_FILE), "assemblies.sampledata")
Utilities.dump_tsv(spadesResults, table_file=spadesResultsSampleData)
print(projectDescriber.SAMPLE_DATA_FILE, "\n", spadesResultsSampleData)
