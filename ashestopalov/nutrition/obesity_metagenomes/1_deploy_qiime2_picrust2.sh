#!/usr/bin/env bash

SAMPLEDATA_MASK="${1}"
if [ -z "${SAMPLEDATA_MASK}" ]
  then
    echo "No arguments supplied"
    exit 1
fi

SAMPLE_SOURCE="$(echo "${SAMPLEDATA_MASK}" | perl -nle 'm/(stool|blood)_chunk_.+/; print $1')"
ROOT_DIR="/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/"
METADATA_TSV="${ROOT_DIR}sample_data/qiime2_meta_data_${SAMPLE_SOURCE}.tsv"

echo Run QIIME2
WORK_DIR="${ROOT_DIR}qiime2"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}" || exit 1
curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/1a_run_qiime2_dada2.sh" \
  -o "run_qiime2_dada2.sh"

force_docker_pull () {
  while true
  do
    if docker pull "$1"
    then
      return
    fi
  done
}

export IMG="qiime2/core:latest"
force_docker_pull "${IMG}"
docker run --rm --net=host \
  -v /data:/data -v /data1:/data1 \
  -e SAMPLEDATA_MASK="${SAMPLEDATA_MASK}" \
  -e METADATA_TSV="${METADATA_TSV}" \
  --workdir="${WORK_DIR}" \
  ${IMG} bash "run_qiime2_dada2.sh"

cd "${ROOT_DIR}" || exit 1



echo Run PICRUSt2
WORK_DIR="${ROOT_DIR}picrust2"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}" || exit 1
curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/1b_run_picrust2.sh" \
  -o "run_picrust2.sh"

export IMG="quay.io/biocontainers/picrust2:2.4.1--py_0"
force_docker_pull "${IMG}"
docker run --rm --net=host \
  -v /data:/data -v /data1:/data1 \
  -e SAMPLEDATA_MASK="${SAMPLEDATA_MASK}" \
  --workdir="${WORK_DIR}" \
  ${IMG} bash "run_picrust2.sh"

cd "${ROOT_DIR}" || exit 1
