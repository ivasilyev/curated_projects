#!/usr/bin/env bash

export REFERENCE_NAME="SILVA"
export REFERENCE_ROOT_DIR="/data/reference/SILVA/"
export REFERENCE_VERSION="138.1"
export SCRIPT_FILE="/tmp/1_fetch_and_import.sh"

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/meta/scripts/silva/1_fetch_and_import.sh" \
    -o "${SCRIPT_FILE}"

REFERENCE_ROOT_DIR="${REFERENCE_ROOT_DIR}" \
REFERENCE_NAME="${REFERENCE_NAME}" \
REFERENCE_VERSION="${REFERENCE_VERSION}" \
bash "${SCRIPT_FILE}"

rm -f "${SCRIPT_FILE}"
