#!/usr/bin/env bash

curl -fsSL \
    "https://gitlab.com/ivasilyev/biological_projects/-/raw/main/ashestopalov/nutrition/mice_fatty_pentylresorcinol/sample_data/main_sampledata.tsv" \
| grep \
    --invert \
    --perl-regexp \
     '(\t14\t|\t60\t)' \
| tee "/tmp/main_sampledata_90.tsv"

# Then push this file into Git
