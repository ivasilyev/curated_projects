#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from meta.templates.ReferenceDescriberTemplate import ReferenceDescriberTemplate


class ReferenceDescriber(ReferenceDescriberTemplate):
    NAME = "hg19"
    VERSION = "2019-07-25-14-42-32"
    ALIAS = "hg19_v2019-07-25-14-42-32"
    DESCRIPTION = "The February 2009 human reference sequence (GRCh37)"
    DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/"
    WEBSITE = "https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg19&redirect=manual&source=genome.ucsc.edu"
    REFDATA = "/data/reference/homo_sapiens/ucsc/hg/hg19/hg19_v2019-07-25-14-42-32/index/hg19_refdata.json"
