#!/usr/bin/env bash

# The SAMPLEDATA_MASK and METADATA_TSV variables are defined externally
mkdir -p "${SAMPLEDATA_MASK}/logs"
cd "${SAMPLEDATA_MASK}" || exit 1

# The output directory must not exist
rm -rf "main_pipeline"

NPROC="$(grep -c '^processor' /proc/cpuinfo)"
Q2DIR="../../qiime2/${SAMPLEDATA_MASK}/"

# From https://github.com/picrust/picrust2/wiki/Workflow
echo Run the PICRUSt2 pipeline
picrust2_pipeline.py --verbose --stratified --coverage --hsp_method mp \
  --processes "${NPROC}" \
  --study_fasta "${Q2DIR}/closed_reference/dna-sequences.fasta" \
  --input "${Q2DIR}/biom/feature-table.biom" \
  --output "main_pipeline" \
  |& tee "logs/picrust2_pipeline.log"

mkdir -p "described/EC_metagenome_out"
echo Convert tables
convert_table.py \
  "main_pipeline/EC_metagenome_out/pred_metagenome_contrib.tsv.gz" \
  --conversion contrib_to_legacy \
  --output "described/EC_metagenome_out/pred_metagenome_contrib.legacy.tsv"

echo Add KEGG ENZYME descriptions
add_descriptions.py -m EC \
  --input "main_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
  --output "described/EC_metagenome_out/pred_metagenome_unstrat_described.tsv"

echo Add KEGG ORTHOLOGY descriptions
add_descriptions.py -m KO \
  --input "main_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
  --output "described/KO_metagenome_out/pred_metagenome_unstrat_described.tsv"

echo Add MetaCyc descriptions
add_descriptions.py -m METACYC \
  --input "main_pipeline/pathways_out/path_abun_unstrat.tsv.gz" \
  --output "described/pathways_out/path_abun_unstrat_described.tsv"

chmod -R 777 "$(pwd)"
cd ..
exit 0
