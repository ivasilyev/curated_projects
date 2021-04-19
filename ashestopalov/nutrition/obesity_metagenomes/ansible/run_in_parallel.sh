#!/usr/bin/env bash
# To be launched from a master node

export ROOT_DIR="/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/"
export QUEUE_FILE="${ROOT_DIR}sample_data/chunks.txt"
export PLAYBOOK_DIR="${HOME}/qiime2_picrust"

mkdir -p "${PLAYBOOK_DIR}"
cd "${PLAYBOOK_DIR}" || exit 1

echo Dry run
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/ansible/dry_run.yml"
ansible-playbook -i "${AWB_HOSTS}" --user "${AWB_UN}" "dry_run.yml"

echo Restore the queue
cp -r "${QUEUE_FILE}.bak" "${QUEUE_FILE}"

echo Main run
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/ansible/main_run.yml"
ansible-playbook -i "${AWB_HOSTS}" --user "${AWB_UN}" "main_run.yml"
