#!/usr/bin/env bash

export ROOT_DIR="/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/"

cd "${ROOT_DIR}" || exit 1
ls "sample_data/chopped" | \
  perl -nle 'm/((stool|blood)[^\.]+)/; print $1' \
  > "sample_data/chunks.txt"
cp "sample_data/chunks.txt" \
  "sample_data/chunks.txt.bak"

while IFS="" read -r p || [ -n "$p" ]
do
  printf '%s\n' "Processing $p"
  bash 1_deploy_qiime2_picrust2.sh "$p"
done < "sample_data/chunks.txt"
