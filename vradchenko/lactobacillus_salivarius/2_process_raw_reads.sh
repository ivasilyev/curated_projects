#!/usr/bin/env bash

# Clean output directory
OUT_DIR="/data1/bio/projects/vradchenko/lactobacillus_salivarius/pga-pe/"
SDATA_DIR="/data1/bio/projects/vradchenko/lactobacillus_salivarius/sample_data/"
SCRIPT="${OUT_DIR}pipeline_handler.py"

# Generate a single sample sampledata
head -n 2 ${SDATA_DIR}raw.sampledata > ${SDATA_DIR}raw.1.sampledata

mkdir -p ${OUT_DIR}
chmod -R 777 ${OUT_DIR}
rm -rf ${OUT_DIR}*

# Download the script file
rm -f ${SCRIPT}
curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/pga-pe-pipeline/pipeline_handler.py" \
  -o ${SCRIPT}

# Run pipeline for a single sample
python3 \
    ${SCRIPT} \
        -i ${SDATA_DIR}raw.1.sampledata \
        -o ${OUT_DIR}1/

# Run pipeline for all samples
python3 \
    ${SCRIPT} \
        -i ${SDATA_DIR}raw.sampledata \
        -f 6 \
        -o ${OUT_DIR}
