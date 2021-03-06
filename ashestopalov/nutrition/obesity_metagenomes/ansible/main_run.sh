#!/usr/bin/env bash

export IMG_QIIME2="qiime2/core:latest"
export IMG_PICRUSt2="quay.io/biocontainers/picrust2:2.4.1--py_0"
export ROOT_DIR="/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/"
export QUEUE_FILE="${ROOT_DIR}sample_data/chunks.txt"

# A 1-60 second random sleep/pause
random_sleep () {
  sleep $((1 + RANDOM % 60))
}

redeploy_script () {
  rm -f "$1"
  while ! [ -s "$1" ]
  do
    curl -fsSL "$2" -o "$1"
  done
}

cd "${ROOT_DIR}" || exit 1

if ! [ -f "${QUEUE_FILE}" ]
then
  print "Not found: ${QUEUE_FILE}"
  exit 1
fi

# Manage logs
LOG_DIR="${ROOT_DIR}pipeline_logs/$(hostname)/"
mkdir -p "${LOG_DIR}"

# The main loop checks if the queue is empty
while true
do
  # Deploy the script
  SCRIPT="${ROOT_DIR}scripts/$(hostname)/deploy_qiime2_picrust2.sh"
  mkdir -p "$(dirname "${SCRIPT}")"
  redeploy_script "${SCRIPT}" "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/1_deploy_qiime2_picrust2.sh"

  random_sleep

  # Verify that the queue file exists and has a size greater than zero
  if ! [ -s "${QUEUE_FILE}" ]
  then
    echo Empty queue: "${QUEUE_FILE}"
    exit 0
  fi

  # Grab the top line of the queue
  ARGS="$(head -n 1 "${QUEUE_FILE}")"
  sed -i '1d' "${QUEUE_FILE}"

  # Run the script
  echo Processing "${ARGS}"
  bash "${SCRIPT}" "${ARGS}" &> "${LOG_DIR}$(hostname)_${ARGS}.log"
done
