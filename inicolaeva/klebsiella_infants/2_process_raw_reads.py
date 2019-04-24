#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# From another console:
export IMG=ivasilyev/spades_cutadapt:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

# Inside of the created container:
cd ${HOME} && \
git clone https://github.com/ivasilyev/curated_projects.git && \
cd ./curated_projects && \
python3
"""

