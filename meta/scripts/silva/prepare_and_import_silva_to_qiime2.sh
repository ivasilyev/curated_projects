#!/usr/bin/env bash

# It is recommended to launch this script into a working QIIME2 environment
# or an instance of the official QIIME2 Docker image, e.g.:
# export IMG="qiime2/core:latest" && docker pull "${IMG}" && docker run --rm --net=host -v /data:/data -it "${IMG}" bash
# cd "/tmp" && curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/meta/scripts/silva/prepare_and_import_silva_to_qiime2.sh"
# bash "prepare_and_import_silva_to_qiime2.sh" |& tee "prepare_and_import_silva_to_qiime2.log"

export REFERENCE_VERSION="138.1"
cd "/data/reference/SILVA/SILVA_v${REFERENCE_VERSION}"

export _BAR="------------------------------------------"
export _ITERATION=0

function xecho () {
    _ITERATION=$((_ITERATION+1))
    printf "\n\n${_BAR}"
    echo "${_ITERATION}. $@"
    echo "${_BAR}"
}


echo "Import the SILVA version ${REFERENCE_VERSION} database assets into QIIME2"
ls -la
# Files required:
#    tax_slv_ssu_${REFERENCE_VERSION}.txt
#    tax_slv_ssu_${REFERENCE_VERSION}.tre
#    taxmap_slv_ssu_ref_nr_${REFERENCE_VERSION}.txt
#    SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_trunc.fasta
#    SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_full_align_trunc.fasta

xecho "Deploy the SILVA DB parser"
git -C "/opt/" clone "https://github.com/mikerobeson/make_SILVA_db.git"
export SDIR="/opt/make_SILVA_db/"

mkdir -p "logs"

xecho "Parse the SILVA taxonomy"
python "${SDIR}parse_silva_taxonomy.py" \
    -s \
    -t "tax_slv_ssu_${REFERENCE_VERSION}.txt" \
    -p "tax_slv_ssu_${REFERENCE_VERSION}.tre" \
    -m "taxmap_slv_ssu_ref_nr_${REFERENCE_VERSION}.txt" \
    -o "SILVA_${REFERENCE_VERSION}_taxonomy.txt" \
|& tee "logs/parse_silva_taxonomy.log"

xecho "Add header to the SILVA taxonomy table"
printf '#OTU ID\ttaxonomy\n' \
| cat - "SILVA_${REFERENCE_VERSION}_taxonomy.txt" \
> "SILVA_${REFERENCE_VERSION}_taxonomy_headed.tsv"

xecho "Remove taxonomy descriptions from FASTA headers, and convert the sequences from RNA to DNA for the unaligned FASTA"
python "${SDIR}convert_rna_to_dna.py" \
    -i "SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_trunc.fasta" \
    -o "SILVA_seqs.fasta" \
|& tee "logs/convert_rna_to_dna_unaligned.log"

xecho "Remove taxonomy descriptions from FASTA headers, and convert the sequences from RNA to DNA for the aligned FASTA"
python "${SDIR}convert_rna_to_dna.py" \
    -g \
    -i "SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_full_align_trunc.fasta" \
    -o "SILVA_align_seqs.fasta" \
|& tee "logs/convert_rna_to_dna_aligned.log"

xecho "Remove sequences containing 5 or more ambiguous bases and/or homopolymers with 8 or more bases"
python "${SDIR}remove_seqs_with_homopolymers.py" \
    -i "SILVA_seqs.fasta" \
    -o "SILVA_seqs_polyfilt.fasta" \
|& tee "logs/remove_seqs_with_homopolymers.log"

xecho "Filter sequences by length, based on taxonomy"
python "${SDIR}filter_seqs_by_length_and_taxonomy.py" \
    -i "SILVA_seqs_polyfilt.fasta" \
    -o "SILVA_seqs_polyfilt_lenfilt.fasta" \
    -t "SILVA_${REFERENCE_VERSION}_taxonomy.txt" \
|& tee "logs/filter_seqs_by_length_and_taxonomy.log"

xecho "Import consensus taxonomy file for the full-length sequences"
qiime tools import \
    --input-path "SILVA_${REFERENCE_VERSION}_taxonomy.txt" \
    --input-format "HeaderlessTSVTaxonomyFormat" \
    --output-path "SILVA-${REFERENCE_VERSION}-full-length-seq-taxonomy.qza" \
    --type 'FeatureData[Taxonomy]' \
|& tee "logs/qiime tools import Taxonomy.log"

xecho "Import FASTA sequences"
qiime tools import \
    --input-path "SILVA_seqs_polyfilt_lenfilt.fasta" \
    --output-path "SILVA-${REFERENCE_VERSION}-SSURef-Full-Seqs.qza" \
    --type 'FeatureData[Sequence]' \
|& tee "logs/qiime tools import Sequence.log"

xecho "Train classifiers for the full-length sequences"
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "SILVA-${REFERENCE_VERSION}-SSURef-Full-Seqs.qza" \
    --i-reference-taxonomy "SILVA-${REFERENCE_VERSION}-full-length-seq-taxonomy.qza" \
    --o-classifier "SILVA-${REFERENCE_VERSION}-SSURef-full-length-classifier.qza" \
|& tee "logs/feature-classifier fit-classifier-naive-bayes.log"

xecho "Import the unrooted tree"
mkdir -p "trees"
qiime tools import \
    --input-path "tax_slv_ssu_${REFERENCE_VERSION}.tre" \
    --output-path "trees/unrooted.qza" \
    --type 'Phylogeny[Unrooted]' \
|& tee "logs/qiime tools import Phylogeny Unrooted.log"

xecho "Root the unrooted tree based on the midpoint rooting method"
qiime phylogeny midpoint-root \
    --verbose \
    --i-tree "trees/unrooted.qza" \
    --o-rooted-tree "trees/rooted.qza" \
|& tee "logs/qiime phylogeny midpoint-root.log"

xecho Export the rooted tree
qiime tools export \
    --input-path "trees/rooted.qza" \
    --output-path "trees" && \
|& tee "logs/qiime tools export rooted trees.log"

mv "trees/tree.nwk" "trees/rooted.nwk"

chmod -R 777 .

echo "The SILVA version ${REFERENCE_VERSION} database assets were imported as QIIME2 artifacts"
