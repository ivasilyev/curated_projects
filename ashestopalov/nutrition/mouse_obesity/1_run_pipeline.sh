#!/usr/bin/env bash

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
    --env SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
    --net=host \
    --interactive \
    --tty "${IMG}" \
        bash -c '
            git pull --quiet;
            python3 ./meta/scripts/qiime2_sample_data.py \
                --extension ".fastq.gz" \
                --input "${RAW_DIR}" \
                --output "${SAMPLEDATA_DIR}"
        '

ROOT_DIR="$(realpath "${ROOT_DIR}")"
SAMPLEDATA_DIR="$(realpath "${SAMPLEDATA_DIR}")"

SCRIPT_DIR="${ROOT_DIR}scripts/"
SAMPLEDATA_CSV="${SAMPLEDATA_DIR}sample_data.csv"
METADATA_TSV="${SAMPLEDATA_DIR}meta_data.tsv"

QIIME2_DIR="${ROOT_DIR}qiime2/"
PICRUST2_DIR="${ROOT_DIR}picrust2/"

echo Run QIIME2
mkdir -p "${QIIME2_DIR}" "${SCRIPT_DIR}"

cd "${SCRIPT_DIR}" || exit 1
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/2_run_qiime2_dada2.sh"

force_docker_pull () {
  while true
  do
    if docker pull "${1}"
    then
      return
    fi
  done
}

cd "${QIIME2_DIR}" || exit 1
export IMG="qiime2/core:latest"
force_docker_pull "${IMG}"
docker run --rm --net=host \
    -v /data:/data \
    -v /data1:/data1 \
    -v /data03:/data03 \
    -v /data04:/data04 \
    -e QIIME2_DIR="${QIIME2_DIR}" \
    -e SAMPLEDATA_CSV="${SAMPLEDATA_CSV}" \
    -e METADATA_TSV="${METADATA_TSV}" \
    --workdir="${QIIME2_DIR}" \
    ${IMG} bash "${SCRIPT_DIR}run_qiime2_dada2.sh"

cd "${ROOT_DIR}" || exit 1



echo Run PICRUSt2
WORK_DIR="${ROOT_DIR}picrust2"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}" || exit 1
curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/3_run_picrust2.sh" \
  -o "${SCRIPT_DIR}run_picrust2.sh"

export IMG="quay.io/biocontainers/picrust2:2.5.0--pyhdfd78af_0"
force_docker_pull "${IMG}"
docker run --rm --net=host \
    -v /data:/data \
    -v /data1:/data1 \
    -v /data03:/data03 \
    -v /data04:/data04 \
    -e QIIME2_DIR="${QIIME2_DIR}" \
    --workdir="${PICRUST2_DIR}" \
    ${IMG} bash "${SCRIPT_DIR}run_picrust2.sh"

cd "${ROOT_DIR}" || exit 1
