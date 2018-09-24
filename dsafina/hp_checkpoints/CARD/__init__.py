#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from meta.scripts.Utilities import Utilities
from dsafina.hp_checkpoints.ProjectDescriber import ProjectDescriber

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3
"""

# Prepare new raw reads for filtering to HG19 DB
raw_reads_dirs = ["/data1/bio/170922/", "/data1/bio/180419_NB501097_0017_AHL55TBGX5/HP/"]
describer = ProjectDescriber()

describer.sampledata = "/data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints_no_hg19.sampledata"

