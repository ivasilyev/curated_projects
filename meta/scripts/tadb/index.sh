#!/usr/bin/env bash

function describe {
    IMG=ivasilyev/curated_projects:latest && \
    docker pull ${IMG} && \
    docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} \
        bash -c '
            git pull --quiet;
            python3 "${HOME}/scripts/curated_projects/meta/scripts/tadb/reference_describer.py" \
                --output "/data/reference/TADB"
        '
}

describe && \
bash "/data/reference/TADB/tadb_v2017.06/index/index.sh" && \
describe &> "/data/reference/TADB/tadb_v2017.06/describe.log"
