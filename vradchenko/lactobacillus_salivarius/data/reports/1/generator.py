#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git pull
LC_ALL=C python3
"""

import os
import jinja2
import pandas as pd
from meta.scripts.Utilities import Utilities
from vradchenko.lactobacillus_salivarius.ProjectDescriber import ProjectDescriber

INDEX_COL_NAME = "sample_name"
SAMPLE_NAMES = ("1sq_FTP", "2sq_FTP", "336g_Nextera", "517_Nextera")
# LOGS_DIR_ROOT = os.path.join(ProjectDescriber.ROOT_DIR, "pga-pe", "log")
# TOOL_NAMES = ("fastqc", "trimmomatic", "cutadapt", "remove_hg", "bowtie2")
TOOL_VERSIONS = dict(fastqc_version="quay.io/biocontainers/fastqc:0.11.8--1",
                     trimmomatic_version="quay.io/biocontainers/trimmomatic:0.39--1",
                     cutadapt_version="quay.io/biocontainers/cutadapt:2.5--py37h516909a_0",
                     bowtie2_version="quay.io/biocontainers/bowtie2:2.3.5--py37he860b03_0",
                     spades_version="quay.io/biocontainers/spades:3.9.1--0")

templates_dir = os.path.join(ProjectDescriber.ROOT_DIR, "reports", "1")
template = jinja2.Template(Utilities.load_string(os.path.join(templates_dir, "template.txt")))

for sample_name in SAMPLE_NAMES:
    # sample_name = SAMPLE_NAMES[0]
    #
    combined_assembly_statistics_df = Utilities.load_tsv(os.path.join(
        ".", ProjectDescriber.OWNER, ProjectDescriber.NAME, "data", "tables", "combined_assembly_statistics.tsv"))
    submission_report_df = Utilities.load_tsv(os.path.join(
        ".", ProjectDescriber.OWNER, ProjectDescriber.NAME, "data", "tables", "ncbi", "submission_report.tsv"))
    #
    submission_combined_df = pd.concat(
        [i.set_index(INDEX_COL_NAME) for i in (combined_assembly_statistics_df, submission_report_df)],
        axis=1, sort=False)
    submission_combined_df.index.names = [INDEX_COL_NAME]
    #
    rendering_dict = submission_combined_df.loc[sample_name, :].to_dict()
    rendering_dict.update(TOOL_VERSIONS)
    #
    out_dir = os.path.join(templates_dir, "out")
    os.makedirs(out_dir, exist_ok=True)
    Utilities.dump_string(template.render(rendering_dict), os.path.join(out_dir, "{}.txt".format(sample_name)))
