#!/usr/bin/env bash

export REFERENCE_VERSION="9.9"
export REFERENCE_DIR="/data03/bio/reference/IGC/igc_v${REFERENCE_VERSION}/"
export REFERENCE_NFASTA="${REFERENCE_DIR}igc_v${REFERENCE_VERSION}.fasta"
export FILE_LIST="${REFERENCE_DIR}file_list.txt."

mkdir -p "${REFERENCE_DIR}"
cd "${REFERENCE_DIR}" || exit 1
rm -rf "${REFERENCE_DIR}"/*

curl -s "https://db.cngb.org/microbiome/genecatalog/genecatalog_human/" \
| grep -oE '(f|ht)tp[s]:\/\/[^ ]+\.gz' \
| sort \
| uniq --unique \
| tee "${FILE_LIST}"

while read URL
    do
    echo "${URL}"
    export FILE="$(
        echo "${URL}" \
        | grep -oE '(f|ht)tp\.[^ ]+\.gz' \
        | perl -nE 'say/([^\/]+\.gz)/'
    )"
    echo Download "${FILE}"
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
    echo Decompress "${FILE}"
    pigz \
        --force \
        --decompress "${REFERENCE_DIR}${FILE}" \
        --processes $(nproc) \
        --recursive \
        --verbose
    done < "${FILE_LIST}"

ln -s \
    "$(
        find \
            "${REFERENCE_DIR}" \
            -regex '.*\.fa$' \
            -type f \
        | head -n 1
    )" \
    "${REFERENCE_NFASTA}"

export IMG="ivasilyev/bwt_filtering_pipeline_worker:latest" && \
docker pull "${IMG}" && \
docker run \
    -e REFERENCE_DIR="${REFERENCE_DIR}" \
    -e REFERENCE_NFASTA="${REFERENCE_NFASTA}" \
    -it \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    "${IMG}" \
    bash -c '
        cd "${REFERENCE_DIR}";
        python3 "${HOME}/scripts/cook_the_reference.py" \
        --input "${REFERENCE_NFASTA}" \
        --output "${REFERENCE_DIR}index";
    '
