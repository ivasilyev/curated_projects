#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from meta.templates.ProjectDescriberTemplate import ProjectDescriberTemplate


class ProjectDescriber(ProjectDescriberTemplate):
    owner = "ndanilova"
    name = "colitis_crohn"
    directory = "/data1/bio/projects/ndanilova/colitis_crohn/"
    # Look "SCFAs_from_KEGG" project for details
    groupdata = "/data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.groupdata"
    sampledata = "/data1/bio/projects/ndanilova/colitis_crohn/colitis_esc_colitis_rem_crohn_esc_crohn_rem_srr.sampledata"
    mask = "no_hg19"


if __name__ == '__main__':
    describer = ProjectDescriber()
