#!/usr/bin/env bash

export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ${IMG} \
    bash -c '
      python3 /home/docker/scripts/cook_the_reference.py \
        --input "/data/reference/MvirDB/mvirdb_v2012.04.28.fasta" \
        --output "/data/reference/MvirDB/mvirdb_v2012.04.28.index";
    '
