#!/usr/bin/env bash

echo Download the Greengenes database assets

export RDIR="Greengenes_v.13.5"
mkdir -p "${RDIR}"
chmod -R 777 "${RDIR}"
cd "${RDIR}"

URL="ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5_otus.tar.gz"
FEXT=$(basename ${URL})
echo "Fetch ${FEXT}"
curl -fsSLO ${URL}
echo "Unpack ${FEXT}"
tar -xvzf "${FEXT}"

echo Import taxonomy file for the full-length sequence
qiime tools import \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path gg_13_5_otus/taxonomy/97_otu_taxonomy.txt \
  --output-path 97_otu-ref-taxonomy-GG.qza \
  --type 'FeatureData[Taxonomy]'

echo Import FASTA sequences
qiime tools import \
  --input-path gg_13_5_otus/rep_set/97_otus.fasta \
  --output-path 97_otus-GG.qza \
  --type 'FeatureData[Sequence]'

echo The Greengenes database assets were imported into QIIME2
