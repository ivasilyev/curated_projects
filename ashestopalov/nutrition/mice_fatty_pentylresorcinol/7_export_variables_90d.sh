#!/usr/bin/env bash

export ROOT_DIR="/data03/bio/projects/ashestopalov/nutrition/mice_fatty_pentylresorcinol_90/"
export RAW_DIR="/data03/bio/rogachev_mice/raw/"

export PIPELINE_DIR="$(realpath "${ROOT_DIR}")/qiime2-picrust2-pipeline/"
export RAW_DIR="$(realpath "${RAW_DIR}")/"
export LOG_DIR="${PIPELINE_DIR}logs/"
export MAIN_SAMPLEDATA_URL="https://gitlab.com/ivasilyev/biological_projects/-/raw/main/ashestopalov/nutrition/mice_fatty_pentylresorcinol/sample_data/main_sampledata_90.tsv"
export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
export SCRIPT_DIR="${PIPELINE_DIR}scripts/"
export PIPELINE_SCRIPT="${SCRIPT_DIR}1_run_pipeline"
export MAIN_SAMPLEDATA_FILE="${SAMPLEDATA_DIR}main_sampledata_90.tsv"

# sudo rm -rf "${ROOT_DIR}"
mkdir \
    --mode 0777 \
    --parents \
    --verbose \
    "${ROOT_DIR}" \
    "${PIPELINE_DIR}" \
    "${LOG_DIR}" \
    "${SCRIPT_DIR}" \
    "${SAMPLEDATA_DIR}"
chmod -R a+rw "${ROOT_DIR}"


export IMG="ivasilyev/curated_projects:latest"
docker pull "${IMG}"
docker run \
    --env MAIN_SAMPLEDATA_URL="${MAIN_SAMPLEDATA_URL}" \
    --env SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
    --net host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    "${IMG}" \
        bash -c '
            git pull --quiet;
            python3 ./meta/sample_data/qiime/split_qiime2_main_sample_data.py \
                --input "${MAIN_SAMPLEDATA_URL}" \
                --output "${SAMPLEDATA_DIR}"
        '

cd "${PIPELINE_DIR}" || exit 1

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/qiime2_picrust2/shell/1_run_pipeline.sh" \
    -o "${PIPELINE_SCRIPT}"

ROOT_DIR="${PIPELINE_DIR}" \
RAW_DIR="${RAW_DIR}" \
SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
bash "${PIPELINE_SCRIPT}" \
|& tee "${LOG_DIR}$(basename "${PIPELINE_SCRIPT}").log"
