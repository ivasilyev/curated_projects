#!/usr/bin/env bash

# It is recommended to launch this script into a working QIIME2 environment
# or an instance of the official QIIME2 Docker image, e.g.:
  export IMG=qiime2/core:latest && \
  docker pull ${IMG} && \
  docker run --rm --net=host -v /data:/data -it ${IMG} bash

echo Import the SILVA database assets into QIIME2
cd /data/reference/SILVA/SILVA_v138
# Files required:
#    tax_slv_ssu_138.txt
#    tax_slv_ssu_138.tre
#    taxmap_slv_ssu_ref_nr_138.txt
#    SILVA_138_SSURef_NR99_tax_silva_trunc.fasta
#    SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta

echo Deploy the SILVA DB parser
git -C "/opt/" clone "https://github.com/mikerobeson/make_SILVA_db.git"
export SDIR="/opt/make_SILVA_db/"

echo Parse the SILVA 138 Taxonomy
python ${SDIR}parse_silva_taxonomy.py -s \
  -t tax_slv_ssu_138.txt \
  -p tax_slv_ssu_138.tre \
  -m taxmap_slv_ssu_ref_nr_138.txt \
  -o SILVA_138_Taxonomy.txt

echo Add header to the SILVA Taxonomy table
printf '#OTU ID\ttaxonomy\n' \
  | cat - SILVA_138_Taxonomy.txt \
  > SILVA_138_Taxonomy_headed.tsv

echo Remove taxonomy descriptions from FASTA headers, and convert the sequences from RNA to DNA for the unaligned FASTA
python ${SDIR}convert_rna_to_dna.py \
  -i SILVA_138_SSURef_NR99_tax_silva_trunc.fasta \
  -o SILVA_seqs.fasta

echo Remove taxonomy descriptions from FASTA headers, and convert the sequences from RNA to DNA for the aligned FASTA
python ${SDIR}convert_rna_to_dna.py -g \
  -i SILVA_138_SSURef_NR99_tax_silva_full_align_trunc.fasta \
  -o SILVA_align_seqs.fasta

echo Remove sequences containing 5 or more ambiguous bases and/or homopolymers with 8 or more bases
python ${SDIR}remove_seqs_with_homopolymers.py \
  -i SILVA_seqs.fasta \
  -o SILVA_seqs_polyfilt.fasta

echo Filter sequences by length, based on taxonomy
python ${SDIR}filter_seqs_by_length_and_taxonomy.py \
    -i SILVA_seqs_polyfilt.fasta \
    -o SILVA_seqs_polyfilt_lenfilt.fasta \
    -t SILVA_138_Taxonomy.txt

echo Import consensus taxonomy file for the full-length sequences
qiime tools import \
  --input-path SILVA_138_Taxonomy.txt  \
  --input-format HeaderlessTSVTaxonomyFormat \
  --output-path Silva-v138-full-length-seq-taxonomy.qza \
  --type 'FeatureData[Taxonomy]'

echo Import FASTA sequences
qiime tools import \
  --input-path SILVA_seqs_polyfilt_lenfilt.fasta \
  --output-path SILVA-138-SSURef-Full-Seqs.qza \
  --type 'FeatureData[Sequence]'

echo Train classifiers for the full-length sequences
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads SILVA-138-SSURef-Full-Seqs.qza \
  --i-reference-taxonomy Silva-v138-full-length-seq-taxonomy.qza \
  --o-classifier SILVA-138-SSURef-full-length-classifier.qza

echo Import the unrooted tree
mkdir -p "trees"
qiime tools import \
  --input-path tax_slv_ssu_138.tre \
  --output-path "trees/unrooted.qza" \
  --type "Phylogeny[Unrooted]"

echo Root the unrooted tree based on the midpoint rooting method
qiime phylogeny midpoint-root --verbose \
  --i-tree "trees/unrooted.qza" \
  --o-rooted-tree "trees/rooted.qza" \
  |& tee "phylogeny midpoint-root.log"

echo Export the rooted tree
qiime tools export \
  --input-path "trees/rooted.qza" \
  --output-path "trees" && \
mv "trees/tree.nwk" "trees/rooted.nwk"

echo The SILVA database assets were imported into QIIME2
