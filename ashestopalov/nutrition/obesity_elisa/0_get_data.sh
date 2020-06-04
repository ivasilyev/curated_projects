#!/usr/bin/env bash

DIRNAME="/data1/bio/projects/ashestopalov/nutrition/obesity_elisa/raw"

mkdir -p ${DIRNAME}
cd ${DIRNAME}


for SRC in "${DIRNAME}"/*
  do
    TGT=$(echo "${SRC}" | iconv -f utf8 -t koi8-r | catdoc -d us-ascii -s koi8-r)
    echo "Rename: '${SRC}' -> '${TGT}'"
    printf "\n"
    mv "${SRC}" "${TGT}"
  done
