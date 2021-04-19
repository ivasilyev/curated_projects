#!/usr/bin/env bash

export IMG_QIIME2="qiime2/core:latest"
export IMG_PICRUSt2="quay.io/biocontainers/picrust2:2.4.1--py_0"
export ROOT_DIR="/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/"
export QUEUE_FILE="${ROOT_DIR}sample_data/chunks.txt"

# Force pull the images
force_docker_pull () {
  while true
  do
    if docker pull "$1"
    then
      return
    fi
  done
}

force_docker_pull "${IMG_QIIME2}"
force_docker_pull "${IMG_PICRUSt2}"

cd "${ROOT_DIR}" || exit 1

if ! [ -f "${QUEUE_FILE}" ]
then
  print "Not found: ${QUEUE_FILE}"
  exit 1
fi

# The main loop checks if the queue is empty
while [ -s "${QUEUE_FILE}" ]
do
  # Grab the top line of the queue
  ARGS="$(head -n 1 "${QUEUE_FILE}")"

  # And remove it from the queue
  sed -i '1d' "${QUEUE_FILE}"

  # Deploy & run the script
  SCRIPT="${ROOT_DIR}scripts/$(hostname)/deploy_qiime2_picrust2.sh"
  mkdir -p "$(dirname "${SCRIPT}")"
  LOG_DIR="${ROOT_DIR}test_pipeline_logs/$(hostname)/"
  mkdir -p "${LOG_DIR}"
  echo "bash ${SCRIPT} ${ARGS}" >> "${LOG_DIR}$(hostname)_${ARGS}.log"
  # A 1-10 second random sleep/pause
  sleep $((1 + RANDOM % 10))
done

echo Empty queue: "${QUEUE_FILE}"
exit 0
