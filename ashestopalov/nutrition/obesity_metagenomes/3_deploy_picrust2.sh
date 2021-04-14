#!/usr/bin/env bash

mkdir -p "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/picrust2"

export IMG=quay.io/biocontainers/picrust2:2.4.1--py_0 && \
docker pull ${IMG} && \
docker run --rm --net=host -it \
  -v /data:/data -v /data1:/data1 \
  -e SAMPLEDATA_MASK="stool_chunk_66" \
  ${IMG} bash
