#!/usr/bin/env bash

mkdir -p "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/qiime2"

export IMG=qiime2/core:latest && \
docker pull ${IMG} && \
docker run --rm --net=host -it \
  -v /data:/data -v /data1:/data1 \
  -e S_CHUNK="${1}" \
  ${IMG} bash -c 'echo ${S_CHUNK}'
