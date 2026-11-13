#!/usr/bin/env bash

# curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/refs/heads/master/ashestopalov/nutrition/elders/1_export_variables.sh" && bash "1_export_variables.sh"

export ROOT_DIR="/data01/bio/projects/ashestopalov/nutrition/elders/"
export RAW_DIR="/data01/bio/projects/ashestopalov/nutrition/elders/raw"

export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
export SCRIPT_DIR="${ROOT_DIR}scripts/"
export PIPELINE_SCRIPT="${SCRIPT_DIR}1_run_pipeline"

export GROUPING_COLUMN_NAME="SubjectID"
export PREV_CONTROL_COLUMN="SubGroup"
export PREV_CONTROL_INDICATOR="Sport"

export PICRUST2_IMG="quay.io/biocontainers/picrust2:2.5.2--pyhdfd78af_0"
export QIIME2_IMG="quay.io/qiime2/core:2023.9"

# sudo rm -rf "${ROOT_DIR}"
mkdir -pv "${ROOT_DIR}" "${SCRIPT_DIR}" "${SAMPLEDATA_DIR}" 
chmod -Rv a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/qiime2_picrust2/shell/1_run_pipeline.sh" \
    -o "${PIPELINE_SCRIPT}"

export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
cd "${SAMPLEDATA_DIR}" || exit 1
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/refs/heads/master/ashestopalov/nutrition/elders/sampledata/qiime2_meta_data.tsv"
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/refs/heads/master/ashestopalov/nutrition/elders/sampledata/qiime2_sample_data.csv"

bash "${PIPELINE_SCRIPT}"

rm -f "$(realpath "${0}")"

