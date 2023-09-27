#!/usr/bin/env bash

export REFERENCE_NAME="SILVA"
export REFERENCE_VERSION="138.1"

export REFERENCE_DIR="/data/reference/${REFERENCE_NAME}/${REFERENCE_NAME}_v${REFERENCE_VERSION}"

echo "Download the ${REFERENCE_NAME} version ${REFERENCE_VERSION} database assets into '${REFERENCE_DIR}'"

export RV="$(echo "${REFERENCE_VERSION}" | sed 's|\.|_|g')"

mkdir -p "${REFERENCE_DIR}"
chmod -R 777 "${REFERENCE_DIR}"
cd "${REFERENCE_DIR}"

for URL in \
    "https://www.arb-silva.de/fileadmin/silva_databases/release_${RV}/Exports/taxonomy/tax_slv_ssu_${REFERENCE_VERSION}.txt.gz" \
    "https://www.arb-silva.de/fileadmin/silva_databases/release_${RV}/Exports/taxonomy/tax_slv_ssu_${REFERENCE_VERSION}.tre.gz" \
    "https://www.arb-silva.de/fileadmin/silva_databases/release_${RV}/Exports/taxonomy/taxmap_slv_ssu_ref_nr_${REFERENCE_VERSION}.txt.gz" \
    "https://www.arb-silva.de/fileadmin/silva_databases/release_${RV}/Exports/SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_trunc.fasta.gz" \
    "https://www.arb-silva.de/fileadmin/silva_databases/release_${RV}/Exports/SILVA_${REFERENCE_VERSION}_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"
    do
        export BASENAME="$(basename ${URL})"
        echo "Fetch '${BASENAME}'"
        export IMG="p3terx/aria2-pro:latest" && \
        docker pull "${IMG}" && \
        docker run \
            --env URL="${URL}" \
            --env REFERENCE_DIR="${REFERENCE_DIR}" \
            --net=host \
            --rm \
            --volume /data:/data \
            --volume /data1:/data1 \
            --volume /data2:/data2 \
            --volume /data03:/data03 \
            --volume /data04:/data04 \
            "${IMG}" bash -c '
                aria2c \
                    --continue \
                    --dir "${REFERENCE_DIR}" \
                    --max-concurrent-downloads $(nproc) \
                    --max-tries 0 \
                    --retry-wait 5 \
                    --split $(nproc) \
                    "${URL}";
                chmod -R a+rw "${REFERENCE_DIR}";
            '
        echo "Unpack ${BASENAME}"
        gunzip "${BASENAME}"
    done

echo "The ${REFERENCE_NAME} version ${REFERENCE_VERSION} database assets were downloaded into '${REFERENCE_DIR}'"
