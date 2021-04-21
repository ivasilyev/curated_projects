#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from meta.templates.ProjectDescriberTemplate import ProjectDescriberTemplate


class ProjectDescriber(ProjectDescriberTemplate):
    OWNER = "ashestopalov"
    NAME = "nutrition/obesity_metagenomes"
    ROOT_DIR = os.path.join("/data1/bio/projects", OWNER, NAME)
    RAW_READS_DIR = os.path.join(ROOT_DIR, "raw_reads")
    SAMPLE_DATA_DIR = os.path.join(ROOT_DIR, "sample_data")
    DATA_DIR = os.path.join(ROOT_DIR, "data")
    MERGED_DATA_DIR = os.path.join(DATA_DIR, "merged_data")
