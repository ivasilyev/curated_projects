#!/usr/bin/env bash

export SRC_DIR="/data1/bio/projects/yutyurin/staphylococcus_aureus_eczema/raw_compressed/"
export ROOT_DIR="/data1/bio/projects/yutyurin/staphylococcus_aureus_eczema/"

mkdir -p "${ROOT_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1

curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/inicolaeva/salmonella_enterica_eclair/1_run_pipeline.sh"

bash "1_run_pipeline.sh"
