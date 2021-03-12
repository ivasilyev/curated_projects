#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from meta.templates.ProjectDescriberTemplate import ProjectDescriberTemplate


class ProjectDescriber(ProjectDescriberTemplate):
    OWNER = "yutyurin"
    NAME = "staphylococcus_aureus"
    ROOT_DIR = os.path.join("/data1/bio/projects", OWNER, NAME)
    RAW_DATA_DIR = os.path.join(ROOT_DIR, "raw")
    MAPPED_DATA_DIR = os.path.join(ROOT_DIR, "mapped")
    DATA_DIGEST_DIR = os.path.join(ROOT_DIR, "digest")
    SAMPLE_DATA_DIR = os.path.join(ROOT_DIR, "sample_data")
    SAMPLE_DATA_FILE = os.path.join(SAMPLE_DATA_DIR, "raw.sampledata")
    GROUP_DATA_FILE = ""
