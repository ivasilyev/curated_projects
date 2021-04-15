#!/usr/bin/env bash

# SAMPLEDATA_MASK="${1}"
SAMPLEDATA_MASK="stool_chunk_1"
METADATA_TSV="/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/sample_data/qiime2_meta_data_$(echo ${SAMPLEDATA_MASK} | grep -E '(stool|blood)').tsv"

echo Run QIIME2
mkdir -p "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/qiime2"
cd "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/qiime2" || exit 1
curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/2_run_qiime2_dada2.sh" \
  -o "run_qiime2_dada2.sh"

export IMG=qiime2/core:latest && \
docker pull ${IMG} && \
docker run --rm --net=host -it \
  -v /data:/data -v /data1:/data1 \
  -e SAMPLEDATA_MASK="${SAMPLEDATA_MASK}" \
  -e METADATA_TSV="${METADATA_TSV}" \
  --workdir="$(pwd)" \
  ${IMG} bash "2_run_qiime2_dada2.sh"



echo Run PICRUSt2
mkdir -p "../picrust2"
cd "../picrust2" || exit 1
curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/4_run_picrust2.sh" \
  -o "run_picrust2.sh"

export IMG=quay.io/biocontainers/picrust2:2.4.1--py_0 && \
docker pull ${IMG} && \
docker run --rm --net=host -it \
  -v /data:/data -v /data1:/data1 \
  -e SAMPLEDATA_MASK="${SAMPLEDATA_MASK}" \
  --workdir="$(pwd)" \
  ${IMG} bash "run_picrust2.sh"
