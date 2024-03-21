#!/usr/bin/env bash

export REFERENCE_ROOT_DIR="/data/reference/SILVA/"
export SCRIPT_FILE="/tmp/1_fetch_and_import.sh"

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/meta/scripts/silva/1_fetch_and_import.sh" \
    -o "${SCRIPT_FILE}"

REFERENCE_ROOT_DIR="${REFERENCE_ROOT_DIR}" \
bash "${SCRIPT_FILE}"

rm -f "${SCRIPT_FILE}"
