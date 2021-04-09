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

echo Run the PICRUSt2 pipeline
mkdir -p "picrust2"
picrust2_pipeline.py --verbose --hsp_method mp \
  --processes ${NPROC} \
  --study_fasta "closed_reference_97/dna-sequences.fasta" \
  --input "biom/feature-table.biom" \
  --output "picrust2"

cd "picrust2"

echo Add KEGG ENZYME descriptions
add_descriptions.py -m EC \
  --input "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
  --output "EC_metagenome_out/pred_metagenome_unstrat_described.tsv"

echo Add KEGG ORTHOLOGY descriptions
add_descriptions.py -m KO \
  --input "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
  --output "KO_metagenome_out/pred_metagenome_unstrat_described.tsv"

echo Add MetaCyc descriptions
add_descriptions.py -m METACYC \
  --input "pathways_out/path_abun_unstrat.tsv.gz" \
  --output "pathways_out/path_abun_unstrat_described.tsv"
