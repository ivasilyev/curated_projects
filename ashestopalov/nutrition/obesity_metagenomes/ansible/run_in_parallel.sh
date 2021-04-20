#!/usr/bin/env bash
# To be launched from a master node

export ROOT_DIR="/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/"
export QUEUE_FILE="${ROOT_DIR}sample_data/chunks.txt"
export PLAYBOOK_DIR="${HOME}/qiime2_picrust2_pipeline"

mkdir -p "${PLAYBOOK_DIR}"
cd "${PLAYBOOK_DIR}" || exit 1

restore_queue() {
  echo Restore the queue
  cp -r "${QUEUE_FILE}.bak" "${QUEUE_FILE}"
}

redeploy_script () {
  rm -f "$1"
  while ! [ -s "$1" ]
  do
    curl -fsSL "$2" -o "$1"
  done
}

echo Dry run
restore_queue
redeploy_script "dry_run.yml" "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/ansible/dry_run.yml"
ansible-playbook -i "${AWB_HOSTS}" --user "${AWB_UN}" "dry_run.yml" -vvvv

echo Main run
restore_queue
redeploy_script "main_run.yml" "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/obesity_metagenomes/ansible/main_run.yml"
rm -f nohup.out
nohup ansible-playbook -i "${AWB_HOSTS}" --user "${AWB_UN}" "main_run.yml" -vvvv
