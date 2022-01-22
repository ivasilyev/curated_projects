#!/usr/bin/env bash

export SRC_DIR="/data1/bio/211022_M04046_0170_000000000-JDB8Y/Conversion_WGS/Eclair_WGS/"
export ROOT_DIR="/data1/bio/projects/inicolaeva/salmonella_enterica_eclair/"

mkdir -p "${ROOT_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1

curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/inicolaeva/salmonella_enterica_eclair/1_run_pipeline.sh"

bash "1_run_pipeline.sh"
