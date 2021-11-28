#!/usr/bin/env bash

function describe {
    IMG=ivasilyev/curated_projects:latest && \
    docker pull ${IMG} && \
    docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} \
        bash -c '
            git pull --quiet;
            python3 "${HOME}/scripts/curated_projects/meta/scripts/mvirdb/reference_describer.py" \
                --output "/data/reference/MvirDB"
        '
}

describe

bash "/data/reference/MvirDB/mvirdb_v2012.04.28/index/index.sh"

describe
