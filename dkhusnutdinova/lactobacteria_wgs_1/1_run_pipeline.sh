#!/usr/bin/env bash

# CARD_VERSION, SRC_DIRS and ROOT_DIR are exported from '0_export_variables.sh'
source "0_export_variables.sh"

export ROOT_DIR="$(realpath -s "${ROOT_DIR}")/"

export RAW_DIR="${ROOT_DIR}raw/"
export SAMPLEDATA_FILE="${ROOT_DIR}sampledata.json"
export PIPELINE_DIR="${ROOT_DIR}pga-pe-pipeline/"


mkdir -p "${RAW_DIR}"
chmod -R a+rw "${RAW_DIR}"
cd "${ROOT_DIR}" || exit 1

for SRC_DIR in ${SRC_DIRS[@]}
    do
        export SRC_DIR="$(realpath -s "${SRC_DIR}")/"
        echo "Unpack reads from ${SRC_DIR}"

        cd "${RAW_DIR}" || exit 1

        find "${SRC_DIR}" \
            -name '*.gz' \
            -type f \
            -print0 \
        | xargs \
            --max-procs "$(nproc)" \
            --null \
            --replace="{}" \
                bash -c '
                    FILE="{}";
                    BASENAME="$(basename "${FILE%.*}")";
                    echo "Decompress ${FILE} to ${RAW_DIR}${BASENAME}";
                    zcat "${FILE}" > "${RAW_DIR}${BASENAME}";
                '
        echo "Deploy symlinks from ${SRC_DIR}"

        find "${SRC_DIR}" \
            -name '*.fastq' \
            -type f \
            -print0 \
        | xargs \
            -0 \
            --max-procs "$(nproc)" \
            -I "{}" \
                bash -c '
                    FILE="{}";
                    BASENAME="$(basename "${FILE}")";
                    echo "Linking ${FILE} to ${RAW_DIR}${BASENAME}";
                    ln -s "${FILE}" "${RAW_DIR}${BASENAME}";
                '
    done


echo "Create sampledata"

export IMG=ivasilyev/curated_projects:latest && \
docker pull "${IMG}" && \
docker run \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
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
    --blast_dir "/data/reference/GenBank" \
    --blast_number 25 \
    --card_json "/data/reference/CARD/card_${CARD_VERSION}/card.json" \
    --hg_dir "/data/reference/homo_sapiens/Ensembl/GRCh38/bowtie2_idx" \
    --input "${SAMPLEDATA_FILE}" \
    --output_dir "${PIPELINE_DIR}" \
    --refdata \
        "/data/reference/CARD/card_${CARD_VERSION}/index/card_${CARD_VERSION}_refdata.json" \
        "/data/reference/MvirDB/mvirdb_v2012.04.28/index/mvirdb_v2012.04.28_refdata.json" \
        "/data/reference/TADB/tadb_v2017.06/index/tadb_v2017.06_refdata.json" \
        "/data/reference/VFDB/vfdb_v2022.01.21/index/vfdb_v2022.01.21_refdata.json" \
    --srst2_dir "/data/reference/SRST2"

rm -f "pipeline_handler.py"



echo "Concatenate coverages"

cd "${PIPELINE_DIR}" || exit 1

export IMG=ivasilyev/curated_projects:latest && \
docker pull "${IMG}" && \/
docker run \
    --env PIPELINE_DIR="${PIPELINE_DIR}" \
    --env CARD_VERSION="${CARD_VERSION}" \
    --interactive \
    --net=host \
    --rm \
    --tty \
    --volume "/data:/data" \
    --volume "/data1:/data1" \
    --volume "/data2:/data2" \
    --volume "/data03:/data03" \
    --volume "/data04:/data04" \
    "${IMG}" \
        bash -c '
            git pull --quiet;
            RGI_DIR="${PIPELINE_DIR}12_rgi/";
            NBEE_DIR="${PIPELINE_DIR}17_nbee_with_annotation/";
            python3 ./mshvydkaya/clostridia_wgs/1_concat_coverage_excel.py \
                --rgi_dir "${RGI_DIR}" \
                --card_version "${CARD_VERSION}" \
                --nbee_dir "${NBEE_DIR}" \
                --output_file "${NBEE_DIR}coverages.xlsx"
        '
