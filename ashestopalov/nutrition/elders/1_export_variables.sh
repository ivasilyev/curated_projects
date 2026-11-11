#!/usr/bin/env bash

#!/usr/bin/env bash

export ROOT_DIR="/data01/bio/projects/ashestopalov/nutrition/elders/"
export RAW_DIR="/data01/bio/projects/ashestopalov/nutrition/elders/raw"

export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
export SCRIPT_DIR="${ROOT_DIR}scripts/"
export PIPELINE_SCRIPT="${SCRIPT_DIR}1_run_pipeline"

# sudo rm -rf "${ROOT_DIR}"
mkdir -p "${ROOT_DIR}" "${SCRIPT_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/qiime2_picrust2/shell/1_run_pipeline.sh" \
    -o "${PIPELINE_SCRIPT}"

export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
cd "${SAMPLEDATA_DIR}" || exit 1
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/refs/heads/master/ashestopalov/nutrition/elders/sampledata/qiime2_meta_data.tsv"
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/refs/heads/master/ashestopalov/nutrition/elders/sampledata/qiime2_sample_data.csv"

bash "${PIPELINE_SCRIPT}"

rm -f "$(realpath "${0}")"

