#!/usr/bin/env bash

export ROOT_DIR="/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/"
export RAW_DIR="/data03/bio/rogachev_human/"



export ROOT_DIR="$(realpath "${ROOT_DIR}")/qiime2-picrust2-pipeline/"
export RAW_DIR="$(realpath "${RAW_DIR}")/"
export LOG_DIR="${ROOT_DIR}logs/"
export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
export SCRIPT_DIR="${ROOT_DIR}scripts/"
export PIPELINE_SCRIPT="${SCRIPT_DIR}1_run_pipeline"

# sudo rm -rf "${ROOT_DIR}"
mkdir -p \
    "${ROOT_DIR}" \
    "${LOG_DIR}" \
    "${SCRIPT_DIR}" \
    "${SAMPLEDATA_DIR}"
chmod -R a+rw "${ROOT_DIR}"

cd "${SAMPLEDATA_DIR}" || exit 1

curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/stool_to_blood_2/sampledata/qiime2_meta_data.tsv"
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/stool_to_blood_2/sampledata/qiime2_sample_data.csv"

cd "${ROOT_DIR}" || exit 1

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/1_run_pipeline.sh" \
    -o "${PIPELINE_SCRIPT}"

ROOT_DIR="${ROOT_DIR}" \
RAW_DIR="${RAW_DIR}" \
SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
bash "${PIPELINE_SCRIPT}" \
|& tee "${LOGS_DIR}$(basename "${PIPELINE_SCRIPT}").log"
