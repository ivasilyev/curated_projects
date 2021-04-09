#!/usr/bin/env bash

export IMG=quay.io/biocontainers/picrust2:2.4.1--py_0 && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 --net=host -it ${IMG} bash

# node5
export SSRC="stool"
# node6
export SSRC="blood"

NPROC="$(grep -c '^processor' /proc/cpuinfo)"
cd "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/qiime2/${SSRC}"

mkdir -p "picrust2"
picrust2_pipeline.py --verbose --hsp_method mp \
  --processes ${NPROC} \
  --study_fasta "closed_reference_97/dna-sequences.fasta" \
  --input "biom/feature-table.biom" \
  --output "picrust2"
