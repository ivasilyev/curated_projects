#!/usr/bin/env bash

echo "Create sampledata"
export IMG=ivasilyev/curated_projects:latest && \
docker pull "${IMG}" && \
docker run \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    --env RAW_DIR="${RAW_DIR}" \
    --env SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
    --net=host \
    --interactive \
    --tty "${IMG}" \
        bash -c '
            git pull --quiet;
            python3 ./meta/scripts/qiime2_sample_data.py \
                --extension ".fastq.gz" \
                --input "${RAW_DIR}" \
                --output "${SAMPLEDATA_DIR}"
        '

ROOT_DIR="$(realpath "${ROOT_DIR}")/"
SAMPLEDATA_DIR="$(realpath "${SAMPLEDATA_DIR}")/"

echo "Working on ${ROOT_DIR}"
SCRIPT_DIR="${ROOT_DIR}scripts/"
SAMPLEDATA_CSV="${SAMPLEDATA_DIR}sample_data.csv"
METADATA_TSV="${SAMPLEDATA_DIR}meta_data.tsv"

QIIME2_DIR="${ROOT_DIR}qiime2/"
PICRUST2_DIR="${ROOT_DIR}picrust2/"

echo "Run QIIME2 in ${QIIME2_DIR}"
mkdir -p "${QIIME2_DIR}" "${PICRUST2_DIR}" "${SCRIPT_DIR}"


force_docker_pull () {
  while true
  do
    if docker pull "${1}"
    then
      return
    fi
  done
}

cd "${SCRIPT_DIR}" || exit 1
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/2_run_qiime2_dada2.sh"
cd "${QIIME2_DIR}" || exit 1
export IMG="qiime2/core:latest"
force_docker_pull "${IMG}"
docker run --rm --net=host \
    -v /data:/data \
    -v /data1:/data1 \
    -v /data03:/data03 \
    -v /data04:/data04 \
    -e QIIME2_DIR="${QIIME2_DIR}" \
    -e SAMPLEDATA_CSV="${SAMPLEDATA_CSV}" \
    -e METADATA_TSV="${METADATA_TSV}" \
    --workdir="${QIIME2_DIR}" \
    ${IMG} bash "${SCRIPT_DIR}2_run_qiime2_dada2.sh"

cd "${ROOT_DIR}" || exit 1



echo "Run PICRUSt2 in ${PICRUST2_DIR}"
cd "${SCRIPT_DIR}" || exit 1
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/3_run_picrust2.sh"
cd "${PICRUST2_DIR}" || exit 1
export IMG="quay.io/biocontainers/picrust2:2.5.0--pyhdfd78af_0"
force_docker_pull "${IMG}"
docker run --rm --net=host \
    -v /data:/data \
    -v /data1:/data1 \
    -v /data03:/data03 \
    -v /data04:/data04 \
    -e QIIME2_DIR="${QIIME2_DIR}" \
    -e PICRUST2_DIR="${PICRUST2_DIR}" \
    --workdir="${PICRUST2_DIR}" \
    ${IMG} bash "${SCRIPT_DIR}3_run_picrust2.sh"

cd "${ROOT_DIR}" || exit 1
