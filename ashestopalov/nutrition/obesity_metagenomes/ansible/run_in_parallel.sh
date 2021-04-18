#!/usr/bin/env bash

export PLAYBOOK_DIR="$HOME/qiime2_picrust"
mkdir -p "${PLAYBOOK_DIR}"
cd "${PLAYBOOK_DIR}" || exit 1

echo Dry run
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/ansible/dry_run.yml"
ansible-playbook -i "${PLAYBOOK_DIR}" --user ${AWB_UN} dry_run.yml



echo Main run
cp -r "${QUEUE_FILE}.bak" "${QUEUE_FILE}"
curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/ansible/main_run.yml"
ansible-playbook -i "${PLAYBOOK_DIR}" --user ${AWB_UN} main_run.yml
