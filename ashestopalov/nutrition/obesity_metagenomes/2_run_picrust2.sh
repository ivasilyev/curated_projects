#!/usr/bin/env bash

mkdir -p "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/picrust2"

export IMG=quay.io/biocontainers/picrust2:2.4.1--py_0 && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 --net=host -it ${IMG} bash

cd "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/picrust2"

# node5
export SSRC="stool"
# node6
export SSRC="blood"

# The output directory must not exist
rm -rf ${SSRC}

NPROC="$(grep -c '^processor' /proc/cpuinfo)"

echo Run the PICRUSt2 pipeline
picrust2_pipeline.py --verbose --stratified  --hsp_method mp \
  --processes "${NPROC}" \
  --study_fasta "closed_reference_97/dna-sequences.fasta" \
  --input "qiime2/${SSRC}/biom/OTU.biom" \
  --output "${SSRC}"


--coverage



chmod -R 777 ${SSRC}
cd "${SSRC}"

metagenome_pipeline.py --strat_out \
  --input "../biom/feature-table.biom" \
  --marker "marker_predicted_and_nsti.tsv.gz" \
  --function "EC_predicted.tsv.gz" \
  --out_dir "EC_metagenome_out"

convert_table.py \
  "EC_metagenome_out/pred_metagenome_contrib.tsv.gz" \
  --conversion contrib_to_legacy \
  --output "EC_metagenome_out/pred_metagenome_contrib.legacy.tsv"

pathway_pipeline.py --verbose \
  --processes "${NPROC}" \
  --input "EC_metagenome_out/pred_metagenome_contrib.tsv.gz" \
  --out_dir "pathways_out"

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
