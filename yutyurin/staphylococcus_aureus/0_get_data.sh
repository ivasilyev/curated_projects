#!/usr/bin/env bash

RAW_DIR="/data1/bio/projects/yutyurin/staphylococcus_aureus/raw/"

mkdir -p ${RAW_DIR}
cp -r "/data2/bio/книэм/"* ${RAW_DIR}
chmod -R 777 /data1/bio/projects/yutyurin

export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git pull && \
jupyter lab --ip=0.0.0.0 --port=61156 --no-browser --NotebookApp.token=TOKEN
