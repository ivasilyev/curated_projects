#!/usr/bin/env bash

export ROOT_DIR="/data01/bio/projects/ashestopalov/nutrition/elders/"
export RAW_DIR="/data01/bio/projects/ashestopalov/nutrition/elders/raw"

export PIPELINE_DIR="$(realpath "${ROOT_DIR}")/qiime2-picrust2-pipeline/"
export RAW_DIR="$(realpath "${RAW_DIR}")/"
export LOG_DIR="${PIPELINE_DIR}logs/"
export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
export SCRIPT_DIR="${PIPELINE_DIR}scripts/"
export PIPELINE_SCRIPT="${SCRIPT_DIR}1_run_pipeline"

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

cd "${SAMPLEDATA_DIR}" || exit 1
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/refs/heads/master/ashestopalov/nutrition/elders/sampledata/qiime2_meta_data.tsv"
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/refs/heads/master/ashestopalov/nutrition/elders/sampledata/qiime2_sample_data.csv"

cd "${PIPELINE_DIR}" || exit 1

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/qiime2_picrust2/shell/1_run_pipeline.sh" \
    -o "${PIPELINE_SCRIPT}"

ROOT_DIR="${PIPELINE_DIR}" \
RAW_DIR="${RAW_DIR}" \
SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
bash "${PIPELINE_SCRIPT}" \
|& tee "${LOG_DIR}$(basename "${PIPELINE_SCRIPT}").log"
