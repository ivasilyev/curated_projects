#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from meta.templates.ProjectDescriberTemplate import ProjectDescriberTemplate
import os


class ProjectDescriber(ProjectDescriberTemplate):
    OWNER = "inicolaeva"
    NAME = "klebsiella_infants"
    ROOT_DIR = os.path.join("/data1/bio/projects", OWNER, NAME)
    SAMPLE_DATA_FILE = os.path.join(ROOT_DIR, "trimmed.sampledata")
    GROUP_DATA_FILE = ""
    RAW_DATA_DIR = os.path.join(ROOT_DIR, "raw")
    DATA_DIGEST_DIR = os.path.join(ROOT_DIR, "digest")
