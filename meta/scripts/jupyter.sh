#!/usr/bin/env bash

# curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/meta/scripts/jupyter.sh"
# bash "jupyter.sh"

force_docker_pull () {
    while true
    do
        if docker pull "${1}"
        then
            return
        fi
    done
}


force_git_pull () {
    while true
    do
        if git pull --quiet
        then
            return
        fi
    done
}


export IMG="ivasilyev/curated_projects:latest" && \
force_docker_pull "${IMG}" && \
docker run \
    --env force_git_pull=force_git_pull \
    --interactive \
    --net=host \
    --rm \
    --tty \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    "${IMG}" bash -c '
        git pull --quiet && \
        jupyter lab \
            --ip=0.0.0.0 \
            --no-browser \
            --NotebookApp.token=TOKEN \
            --port=31522
    '

# Web access: `http://<ip>:31522/?token=TOKEN`
