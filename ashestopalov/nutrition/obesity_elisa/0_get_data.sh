#!/usr/bin/env bash

DIRNAME="/data1/bio/projects/ashestopalov/nutrition/obesity_elisa/raw"

mkdir -p ${DIRNAME}
cd ${DIRNAME} || exit 1


for SRC in "${DIRNAME}"/*
  do
    TGT=$(echo "${SRC}" | iconv -f utf8 -t koi8-r | catdoc -d us-ascii -s koi8-r)
    echo "Rename: '${SRC}' -> '${TGT}'"
    printf "\n"
    mv "${SRC}" "${TGT}"
  done

export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git pull && jupyter lab --ip=0.0.0.0 --port=61156 --no-browser --NotebookApp.token=TOKEN
