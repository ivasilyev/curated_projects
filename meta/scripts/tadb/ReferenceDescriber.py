#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${DOCKER_IMAGE_NAME} python3
"""

from meta.templates.ReferenceDescriberTemplate import ReferenceDescriberTemplate


class ReferenceDescriber(ReferenceDescriberTemplate):
    alias = ""
    name = "TADB"
    description = "An updated database of bacterial type II toxin-antitoxin loci"
    documentation = "https://www.ncbi.nlm.nih.gov/pubmed/?term=29106666"
    url = "http://202.120.12.135/TADB2/index.php"
    refdata = ""


# TODO if db becomes rapidly updating: auto crawl download page: http://202.120.12.135/TADB2/download.html
dl_list = [
    # The experimentally validated Type II TA loci
    "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_T_exp.fas",  # Type II Toxin
    "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_AT_exp.fas",  # Type II Antitoxin (protein)
    "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_RE.fas",  # Type II Regulator (of 3-component type II TA loci)
]
