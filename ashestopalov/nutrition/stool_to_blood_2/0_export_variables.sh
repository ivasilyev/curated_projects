export ROOT_DIR="/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/"
export RAW_DIR="/data03/bio/rogachev_human/"

export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
export SCRIPT_DIR="${ROOT_DIR}scripts/"
export PIPELINE_SCRIPT="${SCRIPT_DIR}1_run_pipeline"

# rm -rf "${ROOT_DIR}"
mkdir -p "${ROOT_DIR}" "${SCRIPT_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/1_run_pipeline.sh" \
    -o "${PIPELINE_SCRIPT}"



ROOT_DIR="${ROOT_DIR}" \
RAW_DIR="${RAW_DIR}" \
SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
bash "${PIPELINE_SCRIPT}"
