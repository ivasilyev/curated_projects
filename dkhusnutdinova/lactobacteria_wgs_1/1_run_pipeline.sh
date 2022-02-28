#!/usr/bin/env bash

export ROOT_DIR="$(realpath -s "${ROOT_DIR}")/"
export SRC_DIR="$(realpath -s "${SRC_DIR}")/"

export RAW_DIR="${ROOT_DIR}raw/"
export SAMPLEDATA_FILE="${ROOT_DIR}sampledata.json"



echo "Deploy symlinks"

mkdir -p "${RAW_DIR}"
chmod -R a+rw "${RAW_DIR}"
cd "${ROOT_DIR}" || exit 1

if [[ "${RAW_DIR}" != "${SRC_DIR}" ]]
    then
    for SRC_READ in "${SRC_DIR}"*.fastq.gz
        do
        BASENAME="$(basename "${SRC_READ}")"
        TGT_READ="${RAW_DIR}${BASENAME}"

        if [ ! -s "${TGT_READ}" ]; then
            echo
            # ln -s "${SRC_READ}" "${TGT_READ}"
        fi
        done
fi



echo "Unpack reads"

cd "${RAW_DIR}"

find "${SRC_DIR}" \
    -name '*.gz' \
    -type f \
    -print0 \
| xargs \
    -0 \
    --max-procs "$(nproc)" \
    -I "{}" \
        bash -c '
            FILE="{}";
            BASENAME="$(basename "${FILE%.*}")";
            echo "Decompress ${FILE} to ${RAW_DIR}${BASENAME}";
            zcat "${FILE}" > "${RAW_DIR}${BASENAME}";
        '

# Sampledata was created manually
