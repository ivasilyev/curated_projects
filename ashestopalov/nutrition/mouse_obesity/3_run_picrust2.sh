#!/usr/bin/env bash

# PICRUST2_DIR and QIIME2_DIR variables are defined externally

export LOGS_DIR="${PICRUST2_DIR}/logs"
mkdir -p "${LOGS_DIR}"

cd "${PICRUST2_DIR}" || exit 1

# The output directory must not exist
export "${PIPELINE_DIR}"="${PICRUST2_DIR}main_pipeline"
rm -rf "${PIPELINE_DIR}"

NPROC="$(grep -c '^processor' /proc/cpuinfo)"

echo Run the PICRUSt2 pipeline
picrust2_pipeline.py --verbose --stratified --coverage --hsp_method mp \
  --processes "${NPROC}" \
  --study_fasta "${QIIME2_DIR}closed_reference/dna-sequences.fasta" \
  --input "${QIIME2_DIR}biom/feature-table.biom" \
  --output "${PIPELINE_DIR}" \
  |& tee "${LOGS_DIR}/picrust2_pipeline.log"

mkdir -p "described/EC_metagenome_out"
echo Convert tables
convert_table.py \
  "${PIPELINE_DIR}/EC_metagenome_out/pred_metagenome_contrib.tsv.gz" \
  --conversion contrib_to_legacy \
  --output "described/EC_metagenome_out/pred_metagenome_contrib.legacy.tsv" \
  |& tee "${LOGS_DIR}/convert_table.log"

echo Add KEGG ENZYME descriptions
add_descriptions.py -m EC \
  --input "${PIPELINE_DIR}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
  --output "described/EC_metagenome_out/pred_metagenome_unstrat_described.tsv"

echo Add KEGG ORTHOLOGY descriptions
add_descriptions.py -m KO \
  --input "${PIPELINE_DIR}/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
  --output "described/KO_metagenome_out/pred_metagenome_unstrat_described.tsv"

echo Add MetaCyc descriptions
add_descriptions.py -m METACYC \
  --input "${PIPELINE_DIR}/pathways_out/path_abun_unstrat.tsv.gz" \
  --output "described/pathways_out/path_abun_unstrat_described.tsv"

chmod -R 777 "$(pwd)"
cd ..
exit 0
