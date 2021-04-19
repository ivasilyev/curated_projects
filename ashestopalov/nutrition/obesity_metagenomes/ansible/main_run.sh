#!/usr/bin/env bash

export IMG_QIIME2="qiime2/core:latest"
export IMG_PICRUSt2="quay.io/biocontainers/picrust2:2.4.1--py_0"
export ROOT_DIR="/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/"
export QUEUE_FILE="${ROOT_DIR}sample_data/chunks.txt"

cd "${ROOT_DIR}" || exit 1

# The main loop checks if the queue is empty
while [ -s "${QUEUE_FILE}" ]
do
  # A 1-10 second random sleep/pause
  sleep $((1 + RANDOM % 10))

  # Grab the top line of the queue
  ARGS="$(head -n 1 "${QUEUE_FILE}")"

  # And remove it from the queue
  sed -i '1d' "${QUEUE_FILE}"

  # Deploy & run the script
  echo Processing "${ARGS}"
  LOG_DIR="${ROOT_DIR}pipeline_logs/$(hostname)/"
  mkdir -p "${LOG_DIR}"
  SCRIPT="/tmp/$(hostname)-deploy_qiime2_picrust2.sh"
  curl -fsSL \
    "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/1_deploy_qiime2_picrust2.sh" \
    -o "${SCRIPT}"
  bash "${SCRIPT}" "${ARGS}" &> "${LOG_DIR}$(hostname)_${ARGS}.log"
done

echo Empty queue: "${QUEUE_FILE}"
exit 0
