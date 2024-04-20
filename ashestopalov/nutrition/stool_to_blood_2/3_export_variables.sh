#!/usr/bin/env bash

export GROUP_NAME="${1}"

if [[ -z "${GROUP_NAME}" ]]
    then
    echo "No group name provided"
    exit 1
    fi



export ROOT_DIR_0="${2}"

if [[ -z "${ROOT_DIR_0}" ]]
    then
    echo "No root directory provided"
    exit 1
    fi



export RAW_DIR="${3}"

if [[ -z "${RAW_DIR}" ]]
    then
    echo "No raw data directory provided"
    exit 1
    fi



export SAMPLEDATA_DIR_0="${4}"

if [[ -z "${SAMPLEDATA_DIR_0}" ]]
    then
    echo "No sample data directory provided"
    exit 1
    fi



export ROOT_DIR="${ROOT_DIR_0}${GROUP_NAME}/"

export ROOT_DIR="$(realpath "${ROOT_DIR}")/qiime2-picrust2-pipeline/"
export RAW_DIR="$(realpath "${RAW_DIR}")/"
export LOG_DIR="${ROOT_DIR}logs/"
export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
export SCRIPT_DIR="${ROOT_DIR}scripts/"
export PIPELINE_SCRIPT="${SCRIPT_DIR}1_run_pipeline"

# sudo rm -rf "${ROOT_DIR}"
mkdir \
    --mode 0777 \
    --parent \
    --verbose \
    "${ROOT_DIR}" \
    "${LOG_DIR}" \
    "${SCRIPT_DIR}" \
    "${SAMPLEDATA_DIR}"

chmod \
    --recursive \
    --verbose \
    a+rw \
    "${ROOT_DIR}"

cd "${SAMPLEDATA_DIR}" || exit 1

ln \
    --verbose \
    --symbolic \
    "${SAMPLEDATA_DIR_0}qiime2_meta_data-${GROUP_NAME}.tsv" \
    "${SAMPLEDATA_DIR}qiime2_meta_data.tsv"

ln \
    --symbolic \
    --verbose \
    "${SAMPLEDATA_DIR_0}qiime2_sample_data-${GROUP_NAME}.csv" \
    "${SAMPLEDATA_DIR}qiime2_sample_data.csv"

cd "${ROOT_DIR}" || exit 1

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/qiime2_picrust2/shell/1_run_pipeline.sh" \
    -o "${PIPELINE_SCRIPT}"

ROOT_DIR="${ROOT_DIR}" \
RAW_DIR="${RAW_DIR}" \
SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
bash "${PIPELINE_SCRIPT}" \
|& tee "${LOGS_DIR}$(basename "${PIPELINE_SCRIPT}").log"
