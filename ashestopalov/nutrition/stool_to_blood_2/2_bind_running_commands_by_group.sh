#!/usr/bin/env bash

export ROOT_DIR="/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/by_group/"
export RAW_DIR="/data03/bio/rogachev_human/"
export SCRIPT_DIR="${ROOT_DIR}qiime2-picrust2-pipeline/scripts/"
export SPLIT_DIR="${ROOT_DIR}qiime2-picrust2-pipeline/sample_data/split/"

# sudo rm -rf "${ROOT_DIR}/qiime2-picrust2-pipeline"

mkdir -p -m 0777 "${SCRIPT_DIR}"

cd "${SCRIPT_DIR}"

curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/stool_to_blood_2/3_export_variables.sh"

while read GROUP_NAME
    do

    echo "${GROUP_NAME}"

    echo "Run Q2P2 pipeline for '${GROUP_NAME}'"

    bash \
        "${SCRIPT_DIR}3_export_variables.sh" \
        "${GROUP_NAME}" \
        "${ROOT_DIR}" \
        "${RAW_DIR}" \
        "${SPLIT_DIR}"

    done < "${SPLIT_DIR}sample_groups.txt"
