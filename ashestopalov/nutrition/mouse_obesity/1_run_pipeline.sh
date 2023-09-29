#!/usr/bin/env bash

function log {
    echo "[$(date '+%d-%m-%Y %H:%M:%S')] $@"
}


force_docker_pull () {
  while true
  do
    if docker pull "${1}"
    then
      return
    fi
  done
}


# Required variables begin
export ROOT_DIR="$(realpath "${ROOT_DIR}")/"
export SAMPLEDATA_DIR="$(realpath "${SAMPLEDATA_DIR}")/"
export SCRIPT_DIR="$(realpath "${SCRIPT_DIR}")/"
# Required variables end

log "Working on ${ROOT_DIR}"
export SAMPLEDATA_CSV="${SAMPLEDATA_DIR}qiime2_sample_data.csv"
export METADATA_TSV="${SAMPLEDATA_DIR}qiime2_meta_data.tsv"

export QIIME2_DIR="${ROOT_DIR}qiime2/"
export QIME2_FEATURES_BIOM="${QIIME2_DIR}bioms/feature-table.biom"
export QIME2_FEATURES_FASTA="${QIIME2_DIR}closed_references/dna-sequences.fasta"
export QIIME2_SCRIPT="${QIIME2_DIR}qiime2.sh"

export REFERENCE_NAME="SILVA"
export REFERENCE_VERSION="138.1"

export REFERENCE_DIR="/data/reference/${REFERENCE_NAME}/${REFERENCE_NAME}_v${REFERENCE_VERSION}"
export TAXA_REFERENCE_FEATURES="/data/reference/${REFERENCE_NAME}/${REFERENCE_NAME}_v${REFERENCE_VERSION}/${REFERENCE_NAME}-v${REFERENCE_VERSION}-full-length-seq-taxonomy.qza"
export TAXA_REFERENCE_CLASSIFIER="/data/reference/${REFERENCE_NAME}/${REFERENCE_NAME}_v${REFERENCE_VERSION}/${REFERENCE_NAME}-${REFERENCE_VERSION}-SSURef-full-length-classifier.qza"
export TAXA_REFERENCE_SEQUENCES="/data/reference/${REFERENCE_NAME}/${REFERENCE_NAME}_v${REFERENCE_VERSION}/${REFERENCE_NAME}-${REFERENCE_VERSION}-SSURef-Full-Seqs.qza"
export TAXA_REFERENCE_HEADER="/data/reference/${REFERENCE_NAME}/${REFERENCE_NAME}_v${REFERENCE_VERSION}/${REFERENCE_NAME}_${REFERENCE_VERSION}_Taxonomy_headed.tsv"

export PICRUST2_DIR="${ROOT_DIR}picrust2/"
export PICRUST2_SCRIPT="${PICRUST2_DIR}picrust2.sh"
export RESULT_DIR="${ROOT_DIR}results/"

export OTU_TABLE="${RESULT_DIR}OTUs_with_taxa.tsv"

log "Check QIIME2 sampledata"
if [ ! -s "${SAMPLEDATA_CSV}" ] && [ ! -s "${METADATA_TSV}" ]
    then
    log "Create sampledata in ${SAMPLEDATA_DIR}"
    export IMG="ivasilyev/curated_projects:latest"
    force_docker_pull "${IMG}"
    docker run \
        --env RAW_DIR="${RAW_DIR}" \
        --env SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
        --net host \
        --rm \
        --volume /data:/data \
        --volume /data1:/data1 \
        --volume /data2:/data2 \
        --volume /data03:/data03 \
        --volume /data04:/data04 \
        "${IMG}" \
            bash -c '
                git pull --quiet;
                python3 ./meta/scripts/qiime2_sample_data.py \
                    --extension ".fastq.gz" \
                    --input "${RAW_DIR}" \
                    --output "${SAMPLEDATA_DIR}"
            '
else
    log "QIIME2 sampledata does exist"
fi

cd "${ROOT_DIR}" || exit 1



mkdir -p "${QIIME2_DIR}" "${PICRUST2_DIR}" "${SCRIPT_DIR}"
curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/2_run_qiime2_dada2.sh" \
    -o "${QIIME2_SCRIPT}"
cd "${QIIME2_DIR}" || exit 1

log "Run QIIME2"
export IMG="qiime2/core:latest"
force_docker_pull "${IMG}"
docker run \
    --env QIIME2_DIR="${QIIME2_DIR}" \
    --env SAMPLEDATA_CSV="${SAMPLEDATA_CSV}" \
    --env METADATA_TSV="${METADATA_TSV}" \
    --env TAXA_REFERENCE_FEATURES="${TAXA_REFERENCE_FEATURES}" \
    --env TAXA_REFERENCE_CLASSIFIER="${TAXA_REFERENCE_CLASSIFIER}" \
    --env TAXA_REFERENCE_SEQUENCES="${TAXA_REFERENCE_SEQUENCES}" \
    --env TAXA_REFERENCE_HEADER="${TAXA_REFERENCE_HEADER}" \
    --net host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    --workdir="${QIIME2_DIR}" \
    "${IMG}" \
    bash "${QIIME2_SCRIPT}"

rm -f "${QIIME2_SCRIPT}"

cd "${ROOT_DIR}" || exit 1



curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/3_run_picrust2.sh" \
    -o "${PICRUST2_SCRIPT}"
cd "${PICRUST2_DIR}" || exit 1

log "Run PICRUSt2"
export IMG="quay.io/biocontainers/picrust2:2.5.0--pyhdfd78af_0"
force_docker_pull "${IMG}"
docker run \
    --env QIME2_FEATURES_BIOM="${QIME2_FEATURES_BIOM}" \
    --env QIME2_FEATURES_FASTA="${QIME2_FEATURES_FASTA}" \
    --env PICRUST2_DIR="${PICRUST2_DIR}" \
    --net host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    --workdir="${PICRUST2_DIR}" \
    "${IMG}" \
    bash "${PICRUST2_SCRIPT}"

rm -f "${PICRUST2_SCRIPT}"

cd "${ROOT_DIR}" || exit 1


log "Copy files"
mkdir -p "${RESULT_DIR}"
find "${ROOT_DIR}" \
    -type f \( \
        -name "path_abun_unstrat_described.tsv" \
        -o -name "pred_metagenome_contrib.legacy.tsv" \
        -o -name "OTUs_with_taxa.tsv" \
    \) -print0 \
    | xargs \
        -0 \
        --max-procs "$(nproc)" \
        -I "{}" \
            bash -c '
                export SRC_FILE="{}";
                echo Copy: "${SRC_FILE}";
                cp -r "${SRC_FILE}" "${RESULT_DIR}$(basename "${SRC_FILE}")"
            '
find "${ROOT_DIR}" \
    -type f \
    -name "pred_metagenome_unstrat_described.tsv" \
    -print0 \
    | xargs \
        -0 \
        --max-procs "$(nproc)" \
        -I "{}" \
            bash -c '
                export SRC_FILE="{}";
                echo Copy: "${SRC_FILE}";
                cp -r "${SRC_FILE}" "${RESULT_DIR}$(basename "$(dirname "${SRC_FILE}")")_$(basename "${SRC_FILE}")"
            '



# The first line of the raw file is '# Constructed from biom file'
sed -i '1d' "${OTU_TABLE}"

log "Concatenate tables"
export IMG="ivasilyev/curated_projects:latest"
force_docker_pull "${IMG}"
docker run \
    --env OTU_TABLE="${OTU_TABLE}" \
    --env TAXA_REFERENCE_HEADER="${TAXA_REFERENCE_HEADER}" \
    --net host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    "${IMG}" \
        bash -c '
            OUT_FILE="${OTU_TABLE%.*}_annotated.tsv";
            echo Concatenate table \"${OTU_TABLE}\" and \"${TAXA_REFERENCE_HEADER}\" into \"${OUT_FILE}\";
            git pull --quiet && \
            python3 ./meta/scripts/concatenate_tables.py \
                --axis 1 \
                --index "#OTU ID" \
                --input \
                    "${TAXA_REFERENCE_HEADER}" \
                    "${OTU_TABLE}" \
                --output "${OUT_FILE}"
        '

log "All pipeline runs ended"
