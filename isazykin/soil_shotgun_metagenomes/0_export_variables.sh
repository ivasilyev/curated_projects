#!/usr/bin/env bash

export ROOT_DIR="/data03/bio/projects/isazykin/soil_shotgun_metagenomes/"
export SRC_DIR="${ROOT_DIR}dnbseq-g400/"

ls "${SRC_DIR}"* | grep -E '(_R1|_R2)' |

mkdir -p "${SRC_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1

# unzip -oq -d "raw_compressed" 'Drop*.zip'

curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/isazykin/soil_shotgun_metagenomes/1_run_pipeline.sh"
bash "1_run_pipeline.sh"
rm -f "1_run_pipeline.sh"
