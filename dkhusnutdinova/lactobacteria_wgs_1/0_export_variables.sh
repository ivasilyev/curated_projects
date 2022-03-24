#!/usr/bin/env bash

export ROOT_DIR="/data04/bio/projects/dkhusnutdinova/lactobacteria_wgs_1/"
export SRC_DIRS=(
  "/data2/bio/Lactobacillus_Dilyara/raw/"
  "/data1/bio/220225_NB501097_0075_AHWH5YBGXF/Conversion_lactobacilli/Lactobacilli/"
)

mkdir -p "${ROOT_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1

curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/dkhusnutdinova/lactobacteria_wgs_1/1_run_pipeline.sh"
bash "1_run_pipeline.sh"
rm -f "1_run_pipeline.sh"

curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/mshvydkaya/clostridia_wgs/2_combine_coverage_excel.sh"
bash "2_combine_coverage_excel.sh"
rm -f "2_combine_coverage_excel.sh"
