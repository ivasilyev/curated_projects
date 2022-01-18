#!/usr/bin/env bash

# SRC_DIR and ROOT_DIR are exported from '0_export_variables.sh'
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

echo "Create sampledata"

export IMG=ivasilyev/curated_projects:latest && \
docker pull "${IMG}" && \
docker run \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --env RAW_DIR="${RAW_DIR}" \
    --env SAMPLEDATA_FILE="${SAMPLEDATA_FILE}" \
    --net=host \
    --interactive \
    --tty "${IMG}" \
        bash -c '
            git pull --quiet;
            python3 ./meta/scripts/sample_data.py \
                --extension ".fastq" \
                --input "${RAW_DIR}" \
                --output "${SAMPLEDATA_FILE}"
        '
# Edit `sampledata.json` if required
# nano "${SAMPLEDATA_FILE}"



echo "Run pipeline"

cd "${ROOT_DIR}" || exit 1
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/pga-pe-pipeline/pipeline_handler.py"

python3 "pipeline_handler.py" \
    --hg_dir "/data/reference/homo_sapiens/Ensembl/GRCh38/bowtie2_idx" \
    --input "${SAMPLEDATA_FILE}" \
    --output_dir "${ROOT_DIR}pga-pe-pipeline" \
    --refdata \
        "/data/reference/MvirDB/mvirdb_v2012.04.28/index/mvirdb_v2012.04.28_refdata.json" \
        "/data/reference/TADB/tadb_v2017.06/index/tadb_v2017.06_refdata.json" \
        "/data/reference/VFDB/vfdb_v2021.09.24/index/vfdb_v2021.09.24_refdata.json"
