#!/usr/bin/env bash

export ROOT_DIR="/data1/bio/projects/yutyurin/staphylococcus_aureus/"

echo Decontaminate submitted assembly, take 1
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash -c \
  '
  git pull;
  python3 $UTILS_DIR/ncbi_contamination_remover.py \
    -i /data1/bio/projects/yutyurin/staphylococcus_aureus/pga-pe/06_plasmid_merger/188staph/188staph_genome.fna \
    -c /data1/bio/projects/yutyurin/staphylococcus_aureus/decontamination/Contamination_1.txt \
    -o /data1/bio/projects/yutyurin/staphylococcus_aureus/decontamination/188staph_genome_decontaminated_1.fna
  '
