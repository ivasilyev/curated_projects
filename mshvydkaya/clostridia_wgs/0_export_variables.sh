#!/usr/bin/env bash

export SRC_DIR="/data1/bio/220225_NB501097_0075_AHWH5YBGXF/Conversion_clostridia/Clostridia/"
export ROOT_DIR="/data03/bio/projects/mshvydkaya/clostridia_wgs_13/"

mkdir -p "${ROOT_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1

curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/inicolaeva/salmonella_enterica_eclair/1_run_pipeline.sh"
bash "1_run_pipeline.sh"
rm -f "1_run_pipeline.sh"
