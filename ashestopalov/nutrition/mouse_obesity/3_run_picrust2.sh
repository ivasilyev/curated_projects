#!/usr/bin/env bash

# PICRUST2_DIR and QIIME2_DIR variables are defined externally

export LOG_DIR="${PICRUST2_DIR}logs/"
export PIPELINE_DIR="${PICRUST2_DIR}main_pipeline/"
export NPROC="$(grep -c '^processor' /proc/cpuinfo)"

mkdir -p "${PICRUST2_DIR}" "${LOG_DIR}"
cd "${PICRUST2_DIR}" || exit 1
echo Ensure that the output directory does not exist
rm -rf "${PIPELINE_DIR}"

echo Run the PICRUSt2 pipeline
picrust2_pipeline.py \
    --coverage \
    --hsp_method mp \
    --input "${QIIME2_DIR}biom/feature-table.biom" \
    --processes "${NPROC}" \
    --study_fasta "${QIIME2_DIR}closed_reference/dna-sequences.fasta" \
    --output "${PIPELINE_DIR}" \
    --stratified \
    --verbose \
    |& tee "${LOG_DIR}picrust2_pipeline.log"

echo Run the PICRUSt2 pathway pipeline
pathway_pipeline.py \
    --input "${PIPELINE_DIR}EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    --intermediate "${PIPELINE_DIR}pathways_out/intermediate" \
    --out_dir "${PIPELINE_DIR}pathways_out" \
    --processes "${NPROC}" \
    --verbose \
    |& tee "${LOG_DIR}picrust2_pathway_pipeline.log"

echo Convert tables
mkdir -p "described/EC_metagenome_out"
convert_table.py \
    "${PIPELINE_DIR}EC_metagenome_out/pred_metagenome_contrib.tsv.gz" \
    --conversion contrib_to_legacy \
    --output "described/EC_metagenome_out/pred_metagenome_contrib.legacy.tsv" \
    |& tee "${LOG_DIR}convert_table.log"

echo Add KEGG ENZYME descriptions
add_descriptions.py \
    --input "${PIPELINE_DIR}EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    --map_type EC \
    --output "described/EC_metagenome_out/pred_metagenome_unstrat_described.tsv"

echo Add KEGG ORTHOLOGY descriptions
add_descriptions.py \
    --input "${PIPELINE_DIR}KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    --map_type KO \
    --output "described/KO_metagenome_out/pred_metagenome_unstrat_described.tsv"

echo Add MetaCyc descriptions
add_descriptions.py -m METACYC \
    --input "${PIPELINE_DIR}pathways_out/path_abun_unstrat.tsv.gz" \
    --output "described/pathways_out/path_abun_unstrat_described.tsv"

chmod -R 777 "$(pwd)"
cd ..
exit 0
