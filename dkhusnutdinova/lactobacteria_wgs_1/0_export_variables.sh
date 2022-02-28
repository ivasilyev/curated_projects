#!/usr/bin/env bash

export SRC_DIR="/data2/bio/Lactobacillus_Dilyara/raw/"
export ROOT_DIR="/data1/bio/projects/dkhusnutdinova/lactobacteria_wgs_1/"

mkdir -p "${ROOT_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1



curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/dkhusnutdinova/lactobacteria_wgs_1/1_run_pipeline.sh"
bash "1_run_pipeline.sh"
rm -f "1_run_pipeline.sh"



echo "Run pipeline"

cd "${ROOT_DIR}" || exit 1
rm -f "pipeline_handler.py"
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/pga-pe-pipeline/pipeline_handler.py"

python3 "pipeline_handler.py" \
    --blast_dir "/data/reference/GenBank" \
    --blast_number 25 \
    --hg_dir "/data/reference/homo_sapiens/Ensembl/GRCh38/bowtie2_idx" \
    --input "${ROOT_DIR}sampledata.json" \
    --output_dir "${ROOT_DIR}pga-pe-pipeline" \
    --refdata \
        "/data/reference/MvirDB/mvirdb_v2012.04.28/index/mvirdb_v2012.04.28_refdata.json" \
        "/data/reference/TADB/tadb_v2017.06/index/tadb_v2017.06_refdata.json" \
        "/data/reference/VFDB/vfdb_v2022.01.21/index/vfdb_v2022.01.21_refdata.json" \
    --srst2_dir "/data/reference/SRST2"

rm -f "pipeline_handler.py"
