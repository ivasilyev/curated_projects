#!/usr/bin/env bash

KTE_DIR="/data1/bio/projects/inicolaeva/klebsiella_infants/kleborate"
rm -rf ${KTE_DIR}
IMG=ivasilyev/kleborate_kaptive:latest && \
docker pull $IMG && \
docker run --rm -v /home:/home --net=host -it $IMG bash -c \
    '
    KTE_DIR="/data1/bio/projects/inicolaeva/klebsiella_infants/kleborate";
    mkdir -p ${KTE_DIR};
    cd ${KTE_DIR};
    Kleborate --all \
        -o ${KTE_DIR}/results.txt \
        -a /data1/bio/projects/inicolaeva/klebsiella_infants/assemblies/*.fna;
    chmod -R 777 ${KTE_DIR}
    '
