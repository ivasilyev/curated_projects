#!/usr/bin/env bash

mkdir -p "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/qiime2"
#  -e S_CHUNK="${1}" \

export IMG=qiime2/core:latest && \
docker pull ${IMG} && \
docker run --rm --net=host -it \
  -v /data:/data -v /data1:/data1 \
  -e SAMPLEDATA_MASK="stool_chunk_66" \
  -e METADATA_TSV="/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/sample_data/qiime2_meta_data_stool.tsv" \
  ${IMG} bash "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/qiime2/2_run_qiime2_dada2.sh"
