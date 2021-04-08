#!/usr/bin/env bash

export IMG=quay.io/biocontainers/picrust2:2.4.1--py_0 && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 --net=host -it ${IMG} bash

# node5
export SSRC="stool"
# node6
export SSRC="blood"

cd "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/qiime2/${SSRC}"


picrust2_pipeline.py --verbose -p 20 --hsp_method mp \
  -s "cr-97/dna-sequences.fasta" \
  -i "feature-table.biom" \
  -o "picrust2_pipeline_out"
