#!/usr/bin/env bash

export ROOT_DIR="/data03/bio/projects/dkhusnutdinova/lactobacteria_wgs_1/"
export CARD_VERSION="v3.2.5"
declare -ax SRC_DIRS=(
    "/data2/bio/Lactobacillus_Dilyara/raw/"
    "/data1/bio/220225_NB501097_0075_AHWH5YBGXF/Conversion_lactobacilli/Lactobacilli/"
    "/data1/bio/220514_M01969_0098_000000000-JWH22/Conversion/lacto/"
    "/data1/bio/221025_M01969_0103_000000000-KFGYJ/Conversion1/Dilyara_lacto/"
)

mkdir -p "${ROOT_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1

curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/dkhusnutdinova/lactobacteria_wgs_1/1_run_pipeline.sh"
bash "1_run_pipeline.sh"
rm -f "1_run_pipeline.sh"
