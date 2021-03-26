#!/usr/bin/env bash

ROOT_DIR="/data1/bio/projects/yutyurin/staphylococcus_aureus/"
PIPELINE_DIR="${ROOT_DIR}pga-pe/"
SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
SCRIPT_EXE="${PIPELINE_DIR}pipeline_handler.py"

echo Clean output directory
mkdir -p ${PIPELINE_DIR}
chmod -R 777 ${ROOT_DIR}
rm -rf ${PIPELINE_DIR}*

echo Download the script file
rm -f ${SCRIPT_EXE}
curl -fsSL "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/pga-pe-pipeline/pipeline_handler.py" \
  -o ${SCRIPT_EXE}

echo Run pipeline for all samples
python3 \
    ${SCRIPT_EXE} \
        -i "${SAMPLEDATA_DIR}raw.sampledata" \
        --hg_dir "/data/reference/homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index" \
        -o ${PIPELINE_DIR}

echo Get relatives for combined assembly
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash -c \
  '
  git pull;
  python3 $UTILS_DIR/blast_nucleotide_sequence.py \
    -i /data1/bio/projects/yutyurin/staphylococcus_aureus/pga-pe/06_plasmid_merger/188staph/188staph_genome.fna \
    -o /data1/bio/projects/yutyurin/staphylococcus_aureus/pga-pe/06_plasmid_merger/188staph
  '
