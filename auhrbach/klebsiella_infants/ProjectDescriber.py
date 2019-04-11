#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from meta.templates.ProjectDescriberTemplate import ProjectDescriberTemplate


class ProjectDescriber(ProjectDescriberTemplate):
    owner = "auhrbach"
    name = "klebsiella_infants"
    directory = "/data1/bio/projects/{}/{}/".format(owner, name)
    groupdata = ()
    sampledata = "{}/main.sampledata".format(directory)
    mask = ""
