#!/usr/bin/env bash

# Required variables start
export REFERENCE_DIR="$(realpath "${REFERENCE_DIR}")/"
export REFERENCE_VERSION="${REFERENCE_VERSION}"
# Required variables end

export RV="$(echo "${REFERENCE_VERSION}" | sed 's|\.|_|g')"
export LOG_DIR="${REFERENCE_DIR}logs/fetch/"
export _START_TIME=$(date +'%s')
export _BAR="------------------------------------------"
export _ITERATION=0


function log() {
    _ITERATION=$((_ITERATION+1))
    printf "\n\n${_BAR}\n\n"
    echo "${_ITERATION}. $@"
    echo "${_BAR}"
}


mkdir \
    --mode 0777 \
    --parents \
    --verbose \
    "${REFERENCE_DIR}"

cd "${REFERENCE_DIR}"

echo "Download the SILVA version ${REFERENCE_VERSION} database assets into '${REFERENCE_DIR}'"

for URL in \
    "https://www.arb-silva.de/fileadmin/silva_databases/release_${RV}/Exports/taxonomy/tax_slv_ssu_${REFERENCE_VERSION}.txt.gz" \
    "https://www.arb-silva.de/fileadmin/silva_databases/release_${RV}/Exports/taxonomy/tax_slv_ssu_${REFERENCE_VERSION}.tre.gz" \
    "https://www.arb-silva.de/fileadmin/silva_databases/release_${RV}/Exports/taxonomy/taxmap_slv_ssu_ref_nr_${REFERENCE_VERSION}.txt.gz" \
    "https://www.arb-silva.de/fileadmin/silva_databases/release_${RV}/Exports/SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_trunc.fasta.gz" \
    "https://www.arb-silva.de/fileadmin/silva_databases/release_${RV}/Exports/SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"
    do
        export BASENAME="$(basename ${URL})"

        export FILE_NAME="${REFERENCE_DIR}${BASENAME}"

        log "Fetch '${BASENAME}'"

        aria2c \
            --continue \
            --dir "${REFERENCE_DIR}" \
            --max-concurrent-downloads $(nproc) \
            --max-tries 0 \
            --retry-wait 5 \
            --split $(nproc) \
            "${URL}"
    done

log "Unpack '${FILE_NAME}'"

find "${SRC_DIR}" \
    -maxdepth 1 \
    -name '*.gz' \
    -type f \
    -print0 \
| xargs \
    -0 \
    --max-procs "$(nproc)" \
    -I "{}" \
        bash -c '
            export FILE_NAME="{}";
            export BASENAME="$(basename "${FILE_NAME%.*}")";
            echo "Decompress ${FILE_NAME} to ${REFERENCE_DIR}${BASENAME}";
            zcat "${FILE_NAME}" > "${REFERENCE_DIR}${BASENAME}";
        '

chmod -R 0777 "${REFERENCE_DIR}";

echo "The SILVA version ${REFERENCE_VERSION} database assets were downloaded into '${REFERENCE_DIR}'"
