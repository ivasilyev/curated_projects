#!/usr/bin/env bash

export IMG="ivasilyev/curated_projects:latest" && \
docker pull "${IMG}" && \
docker run \
    --interactive \
    --net=host \
    --rm \
    --tty \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    "${IMG}" bash

git pull --quiet && \
jupyter lab --ip=0.0.0.0 --port=31522 --no-browser --NotebookApp.token=TOKEN

# Web access: `http://<ip>:31522/?token=TOKEN`
