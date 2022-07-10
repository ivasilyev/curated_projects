#!/usr/bin/env bash

export ROOT_DIR="/data1/bio/projects/ashestopalov/nutrition/mouse_obesity/"
export RAW_DIR="/data03/bio/rogachev_mice/raw/"

export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"

mkdir -p "${ROOT_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1



