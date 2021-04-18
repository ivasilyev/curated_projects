#!/usr/bin/env bash

export PLAYBOOK_DIR="$HOME/qiime2_picrust"
mkdir -p "${PLAYBOOK_DIR}"
cd "${PLAYBOOK_DIR}" || exit 1

echo Dry run
curl -fsSLO ""
ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} dry_run.yml



echo Main run
cp -r "${QUEUE_FILE}.bak" "${QUEUE_FILE}"
curl -fsSLO ""
ansible-playbook -i ${AWB_HOSTS} --user ${AWB_UN} main_run.yml
