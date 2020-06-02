#!/usr/bin/env bash

DIRNAME="/data1/bio/projects/ashestopalov/nutrition/obesity_elisa/raw"

mkdir -p ${DIRNAME}
cd ${DIRNAME}

for FILE in "${DIRNAME}"/*
  do
    HASH=$(openssl sha1 "${FILE}" | grep -o '[^ ]*$')
    EXTENSION=${FILE##*.}
    mv "${FILE}" "${DIRNAME}/${HASH}.${EXTENSION}"
  done
