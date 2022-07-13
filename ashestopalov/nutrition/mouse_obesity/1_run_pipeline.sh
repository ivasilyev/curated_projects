#!/usr/bin/env bash

export ROOT_DIR="$(realpath "${ROOT_DIR}")/"
export SAMPLEDATA_DIR="$(realpath "${SAMPLEDATA_DIR}")/"
export SCRIPT_DIR="$(realpath "${SCRIPT_DIR}")/"

echo "Working on ${ROOT_DIR}"
export SAMPLEDATA_CSV="${SAMPLEDATA_DIR}sample_data.csv"
export METADATA_TSV="${SAMPLEDATA_DIR}meta_data.tsv"

export QIIME2_DIR="${ROOT_DIR}qiime2/"
export PICRUST2_DIR="${ROOT_DIR}picrust2/"


force_docker_pull () {
  while true
  do
    if docker pull "${1}"
    then
      return
    fi
  done
}


echo "Create sampledata in ${SAMPLEDATA_DIR}"
export IMG="ivasilyev/curated_projects:latest"
force_docker_pull "${IMG}"
docker run \
    --env RAW_DIR="${RAW_DIR}" \
    --env SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
    --net=host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    "${IMG}" \
        bash -c '
            git pull --quiet;
            python3 ./meta/scripts/qiime2_sample_data.py \
                --extension ".fastq.gz" \
                --input "${RAW_DIR}" \
                --output "${SAMPLEDATA_DIR}"
        '

cd "${ROOT_DIR}" || exit 1



mkdir -p "${QIIME2_DIR}" "${PICRUST2_DIR}" "${SCRIPT_DIR}"
export SCRIPT_FILE="${SCRIPT_DIR}2_run_qiime2_dada2.sh"
curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/2_run_qiime2_dada2.sh" \
    -o "${SCRIPT_FILE}"
cd "${QIIME2_DIR}" || exit 1

export IMG="qiime2/core:latest"
force_docker_pull "${IMG}"
docker run \
    --env QIIME2_DIR="${QIIME2_DIR}" \
    --env SAMPLEDATA_CSV="${SAMPLEDATA_CSV}" \
    --env METADATA_TSV="${METADATA_TSV}" \
    --net=host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    --workdir="${QIIME2_DIR}" \
    "${IMG}" bash "${SCRIPT_FILE}"

cd "${ROOT_DIR}" || exit 1



export SCRIPT_FILE="${SCRIPT_DIR}3_run_picrust2.sh"
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/3_run_picrust2.sh" \
    -o "${SCRIPT_FILE}"

export IMG="quay.io/biocontainers/picrust2:2.5.0--pyhdfd78af_0"
force_docker_pull "${IMG}"
docker run \
    --env QIIME2_DIR="${QIIME2_DIR}" \
    --env PICRUST2_DIR="${PICRUST2_DIR}" \
    --net=host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    --workdir="${PICRUST2_DIR}" \
    "${IMG}" bash "${SCRIPT_FILE}"

cd "${ROOT_DIR}" || exit 1
