#!/usr/bin/env bash

export ROOT_DIR="$(realpath "${ROOT_DIR}")/"
export SAMPLEDATA_DIR="$(realpath "${SAMPLEDATA_DIR}")/"
export SCRIPT_DIR="$(realpath "${SCRIPT_DIR}")/"

echo "Working on ${ROOT_DIR}"
export SAMPLEDATA_CSV="${SAMPLEDATA_DIR}sample_data.csv"
export METADATA_TSV="${SAMPLEDATA_DIR}meta_data.tsv"

export QIIME2_DIR="${ROOT_DIR}qiime2/"
export PICRUST2_DIR="${ROOT_DIR}picrust2/"

export QIIME2_SCRIPT="${QIIME2_DIR}qiime2.sh"
export PICRUST2_SCRIPT="${PICRUST2_DIR}picrust2.sh"

export RESULT_DIR="${ROOT_DIR}results/"
export OTU_TABLE="${RESULT_DIR}OTUs.tsv"

force_docker_pull () {
  while true
  do
    if docker pull "${1}"
    then
      return
    fi
  done
}


echo "Create sampledata in ${SAMPLEDATA_DIR}"
export IMG="ivasilyev/curated_projects:latest"
force_docker_pull "${IMG}"
docker run \
    --env RAW_DIR="${RAW_DIR}" \
    --env SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
    --net=host \
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

cd "${ROOT_DIR}" || exit 1



mkdir -p "${QIIME2_DIR}" "${PICRUST2_DIR}" "${SCRIPT_DIR}"
curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/2_run_qiime2_dada2.sh" \
    -o "${QIIME2_SCRIPT}"
cd "${QIIME2_DIR}" || exit 1

echo "Run QIIME2"
export IMG="qiime2/core:latest"
force_docker_pull "${IMG}"
docker run \
    --env QIIME2_DIR="${QIIME2_DIR}" \
    --env SAMPLEDATA_CSV="${SAMPLEDATA_CSV}" \
    --env METADATA_TSV="${METADATA_TSV}" \
    --net=host \
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

echo "Run PICRUSt2"
export IMG="quay.io/biocontainers/picrust2:2.5.0--pyhdfd78af_0"
force_docker_pull "${IMG}"
docker run \
    --env QIIME2_DIR="${QIIME2_DIR}" \
    --env PICRUST2_DIR="${PICRUST2_DIR}" \
    --env PICRUST2_SCRIPT="${PICRUST2_SCRIPT}" \
    --net=host \
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


echo "Copy files"
mkdir -p "${RESULT_DIR}"
find "${ROOT_DIR}" \
    -type f \( \
        -name "path_abun_unstrat_described.tsv" \
        -o -name "pred_metagenome_contrib.legacy.tsv" \
        -o -name "OTUs.tsv" \
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

echo "Concatenate tables"
export IMG="ivasilyev/curated_projects:latest"
force_docker_pull "${IMG}"
docker run \
    --env OTU_TABLE="${OTU_TABLE}" \
    --net=host \
    --rm \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    "${IMG}" \
        bash -c '
            git pull --quiet;
            python3 ./meta/scripts/concatenate_tables.py \
                --axis 1 \
                --index "#OTU ID" \
                --input "${OTU_TABLE}" "/data/reference/SILVA/SILVA_v138/SILVA_138_Taxonomy_headed.tsv" / \
                --output "${OTU_TABLE%.*}_annotated.tsv"
        '

echo "All pipeline runs ended"
