#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from meta.templates.ProjectDescriberTemplate import ProjectDescriberTemplate
import os


class ProjectDescriber(ProjectDescriberTemplate):
    owner = "auhrbach"
    name = "klebsiella_infants"
    directory = os.path.join("/data1/bio/projects", owner, name)
    groupdata = ()
    sampledata = os.path.join(directory, "main.sampledata")
