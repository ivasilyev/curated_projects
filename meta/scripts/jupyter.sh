#!/usr/bin/env bash

export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git pull && jupyter lab --ip=0.0.0.0 --port=31522 --no-browser --NotebookApp.token=TOKEN

# Web access: http:<ip>:31522/?token=TOKEN
