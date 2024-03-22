#!/usr/bin/env bash

# It is recommended to launch this script into a working QIIME2 environment
# or an instance of the official QIIME2 Docker image, e.g.:
# export IMG="qiime2/core:latest" && docker pull "${IMG}" && docker run --rm --net=host -v /data:/data -it "${IMG}" bash
# cd "/tmp" && curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/meta/scripts/silva/prepare_and_import_silva_to_qiime2.sh"
# bash "prepare_and_import_silva_to_qiime2.sh" |& tee "prepare_and_import_silva_to_qiime2.log"

# Required variables start
export REFERENCE_DIR="$(realpath "${REFERENCE_DIR}")/"
export REFERENCE_VERSION="${REFERENCE_VERSION}"
# Required variables end

export LOG_DIR="${REFERENCE_DIR}logs/import/"
export _LOG_BAR="=========================================="
export _LOG_COUNTER=1


function log() {
    printf "\n${_LOG_BAR}\n\n[$(date '+%d-%m-%Y %H:%M:%S.%N')][Import][OP#${_LOG_COUNTER}] $@\n\n${_LOG_BAR}\n\n"
    _LOG_COUNTER=$((_LOG_COUNTER + 1))
}


mkdir \
    --mode 0777 \
    --parents \
    --verbose \
    "${LOG_DIR}" \
    "/opt/"

cd "/opt/"

log "Deploy the SILVA DB parser"

git \
    -C "/opt/" \
    clone \
    --verbose \
    "https://github.com/mikerobeson/make_SILVA_db.git"

export PARSER_DIR="/opt/make_SILVA_db/"

cd "${REFERENCE_DIR}"

log "Import the SILVA version ${REFERENCE_VERSION} database assets as QIIME2 into '${REFERENCE_DIR}'"

log "View the directory '${REFERENCE_DIR}' contents"

ls -la
# Files required:
#    tax_slv_ssu_${REFERENCE_VERSION}.txt
#    tax_slv_ssu_${REFERENCE_VERSION}.tre
#    taxmap_slv_ssu_ref_nr_${REFERENCE_VERSION}.txt
#    SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_trunc.fasta
#    SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_full_align_trunc.fasta


log "Parse the SILVA taxonomy"

python "${PARSER_DIR}parse_silva_taxonomy.py" \
    -s \
    -t "${REFERENCE_DIR}tax_slv_ssu_${REFERENCE_VERSION}.txt" \
    -p "${REFERENCE_DIR}tax_slv_ssu_${REFERENCE_VERSION}.tre" \
    -m "${REFERENCE_DIR}taxmap_slv_ssu_ref_nr_${REFERENCE_VERSION}.txt" \
    -o "${REFERENCE_DIR}SILVA_${REFERENCE_VERSION}_taxonomy.txt" \
|& tee "${LOG_DIR}parse_silva_taxonomy.log"



log "Add header to the SILVA taxonomy table"

printf '#OTU ID\ttaxonomy\n' \
| cat - "${REFERENCE_DIR}SILVA_${REFERENCE_VERSION}_taxonomy.txt" \
> "SILVA_${REFERENCE_VERSION}_taxonomy_headed.tsv"



log "Remove taxonomy descriptions from FASTA headers, and convert the sequences from RNA to DNA for the unaligned FASTA"

python "${PARSER_DIR}convert_rna_to_dna.py" \
    -i "${REFERENCE_DIR}SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_trunc.fasta" \
    -o "${REFERENCE_DIR}SILVA_seqs.fasta" \
|& tee "${LOG_DIR}convert_rna_to_dna_unaligned.log"



log "Remove taxonomy descriptions from FASTA headers, and convert the sequences from RNA to DNA for the aligned FASTA"

python "${PARSER_DIR}convert_rna_to_dna.py" \
    -g \
    -i "${REFERENCE_DIR}SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_full_align_trunc.fasta" \
    -o "${REFERENCE_DIR}SILVA_align_seqs.fasta" \
|& tee "${LOG_DIR}convert_rna_to_dna_aligned.log"



log "Remove sequences containing 5 or more ambiguous bases and/or homopolymers with 8 or more bases"

python "${PARSER_DIR}remove_seqs_with_homopolymers.py" \
    -i "${REFERENCE_DIR}SILVA_seqs.fasta" \
    -o "${REFERENCE_DIR}SILVA_seqs_polyfilt.fasta" \
|& tee "${LOG_DIR}remove_seqs_with_homopolymers.log"



log "Filter sequences by length, based on taxonomy"

python "${PARSER_DIR}filter_seqs_by_length_and_taxonomy.py" \
    -i "${REFERENCE_DIR}SILVA_seqs_polyfilt.fasta" \
    -o "${REFERENCE_DIR}SILVA_seqs_polyfilt_lenfilt.fasta" \
    -t "${REFERENCE_DIR}SILVA_${REFERENCE_VERSION}_taxonomy.txt" \
|& tee "${LOG_DIR}filter_seqs_by_length_and_taxonomy.log"



log "Import consensus taxonomy file for the full-length sequences"

qiime tools import \
    --input-path "${REFERENCE_DIR}SILVA_${REFERENCE_VERSION}_taxonomy.txt" \
    --input-format "HeaderlessTSVTaxonomyFormat" \
    --output-path "${REFERENCE_DIR}SILVA-${REFERENCE_VERSION}-full-length-seq-taxonomy.qza" \
    --type 'FeatureData[Taxonomy]' \
|& tee "${LOG_DIR}qiime tools import Taxonomy.log"



log "Import FASTA sequences"

qiime tools import \
    --input-path "${REFERENCE_DIR}SILVA_seqs_polyfilt_lenfilt.fasta" \
    --output-path "${REFERENCE_DIR}SILVA-${REFERENCE_VERSION}-SSURef-Full-Seqs.qza" \
    --type 'FeatureData[Sequence]' \
|& tee "${LOG_DIR}qiime tools import Sequence.log"



log "Train classifiers for the full-length sequences"

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "${REFERENCE_DIR}SILVA-${REFERENCE_VERSION}-SSURef-Full-Seqs.qza" \
    --i-reference-taxonomy "${REFERENCE_DIR}SILVA-${REFERENCE_VERSION}-full-length-seq-taxonomy.qza" \
    --o-classifier "${REFERENCE_DIR}SILVA-${REFERENCE_VERSION}-SSURef-full-length-classifier.qza" \
|& tee "${LOG_DIR}feature-classifier fit-classifier-naive-bayes.log"



log "Import the unrooted tree"

mkdir -p "${REFERENCE_DIR}trees"

qiime tools import \
    --input-path "${REFERENCE_DIR}tax_slv_ssu_${REFERENCE_VERSION}.tre" \
    --output-path "${REFERENCE_DIR}trees/unrooted.qza" \
    --type 'Phylogeny[Unrooted]' \
|& tee "${LOG_DIR}qiime tools import Phylogeny Unrooted.log"



log "Root the unrooted tree based on the midpoint rooting method"

qiime phylogeny midpoint-root \
    --verbose \
    --i-tree "${REFERENCE_DIR}trees/unrooted.qza" \
    --o-rooted-tree "${REFERENCE_DIR}trees/rooted.qza" \
|& tee "${LOG_DIR}qiime phylogeny midpoint-root.log"



log "Export the unrooted trees"

qiime tools export \
    --input-path "${REFERENCE_DIR}trees/unrooted.qza" \
    --output-path "${REFERENCE_DIR}trees" \
|& tee "${LOG_DIR}qiime tools export unrooted trees.log"

mv \
    "${REFERENCE_DIR}trees/tree.nwk" \
    "${REFERENCE_DIR}trees/unrooted.nwk"



log "Export the rooted trees"

qiime tools export \
    --input-path "${REFERENCE_DIR}trees/rooted.qza" \
    --output-path "${REFERENCE_DIR}trees" \
|& tee "${LOG_DIR}qiime tools export rooted trees.log"

mv "${REFERENCE_DIR}trees/tree.nwk" "${REFERENCE_DIR}trees/rooted.nwk"

chmod \
    --verbose \
    --recursive \
    0777 \
    "${REFERENCE_DIR}"

echo "The SILVA version ${REFERENCE_VERSION} database assets were imported as QIIME2 artifacts into '${REFERENCE_DIR}'"
