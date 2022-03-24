#!/usr/bin/env bash

export PIPELINE_DIR="${ROOT_DIR}pga-pe-pipeline"
cd "${PIPELINE_DIR}" || exit 1

export IMG=ivasilyev/curated_projects:latest && \
docker pull "${IMG}" && \
docker run \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    --env PIPELINE_DIR="${PIPELINE_DIR}" \
    --net=host \
    --interactive \
    --tty "${IMG}" \
        bash -c '
            git pull --quiet;
            python3 ./mshvydkaya/clostridia_wgs/1_concat_coverage_excel.py \
                --rgi_dir "${PIPELINE_DIR}/12_rgi" \
                --card_version "v.3.1.4" \
                --nbee_dir "${PIPELINE_DIR}/17_nbee_with_annotation" \
                --output_file "${PIPELINE_DIR}/coverages.xlsx"
        '
