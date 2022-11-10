#!/usr/bin/env bash

export ROOT_DIR="/data03/bio/projects/ashestopalov/nutrition/mouse_obesity/"
export RAW_DIR="/data03/bio/rogachev_mice/raw/"

export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
export SCRIPT_DIR="${ROOT_DIR}scripts/"
export PIPELINE_SCRIPT="${SCRIPT_DIR}1_run_pipeline"

# rm -rf "${ROOT_DIR}"
mkdir -p "${ROOT_DIR}" "${SCRIPT_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/1_run_pipeline.sh" \
    -o "${PIPELINE_SCRIPT}"

bash "${PIPELINE_SCRIPT}"
