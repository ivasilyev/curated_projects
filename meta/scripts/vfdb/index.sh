#!/usr/bin/env bash

function describe {
    IMG=ivasilyev/curated_projects:latest && \
    docker pull ${IMG} && \
    docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} \
        bash -c '
            git pull --quiet;
            python3 "${HOME}/scripts/curated_projects/meta/scripts/vfdb/reference_describer.py" \
                --output "/data/reference/VFDB"
        '
}

LOG="/data/reference/VFDB/describe.log"
describe | tee -a "${LOG}" && \
bash "/data/reference/VFDB/vfdb_v2022.01.21/index/index.sh" | tee -a "${LOG}" && \
describe | tee -a "${LOG}"
