#!/usr/bin/env bash

ROOT_DIR="/data1/bio/projects/yutyurin/staphylococcus_aureus/"
PIPELINE_DIR="${ROOT_DIR}pga-pe/"
SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
SCRIPT_EXE="${PIPELINE_DIR}pipeline_handler.py"

echo Clean output directory
mkdir -p ${PIPELINE_DIR}
chmod -R 777 ${ROOT_DIR}
rm -rf ${PIPELINE_DIR}*

echo Download the script file
rm -f ${SCRIPT_EXE}
curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/pga-pe-pipeline/pipeline_handler.py" \
  -o ${SCRIPT_EXE}

echo Run pipeline for all samples
python3 \
    ${SCRIPT_EXE} \
        -i "${SAMPLEDATA_DIR}raw.sampledata" \
        --hg_dir "/data/reference/homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index" \
        -o ${PIPELINE_DIR}


