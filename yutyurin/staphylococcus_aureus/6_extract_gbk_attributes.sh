#!/usr/bin/env bash

OUT="../isolation_sources.tsv"
printf "geninfo_id\tgenbank_id\tisolation_source\thost\n" > "${OUT}"
for i in *.gbk
    do
        printf "$(basename "${i%.*}")\t$(grep -oP '(?<=VERSION).+$' "${i}" | xargs)\t$(grep -oP '(?<=\/isolation_source=\")[^\"]+(?=\")' "${i}")\t$(grep -oP '(?<=\/host=\")[^\"]+(?=\")' "${i}")\n" >> "${OUT}"
    done
