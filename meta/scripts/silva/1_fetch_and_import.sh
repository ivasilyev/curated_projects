#!/usr/bin/env bash

# Required variables start
export REFERENCE_NAME="${REFERENCE_NAME}"
export REFERENCE_ROOT_DIR="$(realpath "${REFERENCE_ROOT_DIR}")/"
export REFERENCE_VERSION="${REFERENCE_VERSION}"
# Required variables end

export LOG_DIR="${REFERENCE_DIR}logs/pipeline/"
export _LOG_BAR="=========================================="
export _LOG_COUNTER=1


function log() {
    printf "\n${_LOG_BAR}\n\n[$(date '+%d-%m-%Y %H:%M:%S.%N')][Pipeline][OP#${_LOG_COUNTER}] $@\n\n${_LOG_BAR}\n\n"
    _LOG_COUNTER=$((_LOG_COUNTER + 1))
}


export SCRIPT_FILE="/tmp/script.sh"

export REFERENCE_DIR="${REFERENCE_ROOT_DIR}/${REFERENCE_NAME}_v${REFERENCE_VERSION}"

# sudo rm -rf "${REFERENCE_DIR}"


force_docker_pull () {
  while true
  do
    if docker pull "${1}"
    then
      return
    fi
  done
}


get_latest_quay_tag() {
    echo "$(
        curl -fsSL "https://quay.io/api/v1/repository/${1}/${2}" \
        | grep \
            --only-matching \
            --perl-regexp \
            '(?<=\"tags\": {\")[^\"]*(?=\":)'
    )"
}


compose_quay_img() {
    echo "quay.io/${1}/${2}:$(get_latest_quay_tag ${1} ${2})"
}


mkdir \
    --mode 0777 \
    --parents \
    --verbose \
    "${REFERENCE_DIR}"

chmod -R 0777 "${REFERENCE_DIR}"

cd "${REFERENCE_DIR}"



rm -rf "${SCRIPT_FILE}"

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/meta/scripts/silva/2_fetch.sh" \
    -o "${SCRIPT_FILE}"

export IMG="p3terx/aria2-pro:latest"

docker pull "${IMG}"

docker run \
    --env URL="${URL}" \
    --env REFERENCE_DIR="${REFERENCE_DIR}" \
    --env REFERENCE_VERSION="${REFERENCE_VERSION}" \
    --env SCRIPT_FILE="${SCRIPT_FILE}" \
    --net=host \
    --rm \
    --volume "${REFERENCE_DIR}":"${REFERENCE_DIR}" \
    --volume "${SCRIPT_FILE}":"${SCRIPT_FILE}" \
    "${IMG}" bash "${SCRIPT_FILE}"

rm -f "${SCRIPT_FILE}"



curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/meta/scripts/silva/3_import.sh" \
    -o "${SCRIPT_FILE}"

export IMG="$(compose_quay_img qiime2 core)"

force_docker_pull "${IMG}"

docker run \
    --env URL="${URL}" \
    --env REFERENCE_DIR="${REFERENCE_DIR}" \
    --env REFERENCE_VERSION="${REFERENCE_VERSION}" \
    --env SCRIPT_FILE="${SCRIPT_FILE}" \
    --net=host \
    --rm \
    --volume "${REFERENCE_DIR}":"${REFERENCE_DIR}" \
    --volume "${SCRIPT_FILE}":"${SCRIPT_FILE}" \
    "${IMG}" bash "${SCRIPT_FILE}"

rm -f "${SCRIPT_FILE}"

chmod \
    --verbose \
    --recursive \
    0777 \
    "${REFERENCE_DIR}"
