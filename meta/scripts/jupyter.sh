#!/usr/bin/env bash

# curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/meta/scripts/jupyter.sh"
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

## Pre-setup for SSH-based hairpin control
#cat <<EOF | tee -a "${HOME}/.ssh/config"
#
#Host localhost
#    HostName localhost
#    Port 22
#    User ${USER}
#EOF
#ssh-copy-id -i ~/.ssh/id_rsa.pub -o 'StrictHostKeyChecking=no' localhost
#ssh localhost "echo \$(whoami)"

export IMG="ivasilyev/curated_projects:latest" && \
force_docker_pull "${IMG}" && \
docker run \
    --env force_git_pull=force_git_pull \
    --env PORT=31522 \
    --env TOKEN=TOKEN \
    --env IP_ADDRESS="$(ip route | grep default | awk '{ print $3 }')" \
    --interactive \
    --net=host \
    --rm \
    --tty \
    --volume /data:/data \
    --volume /data1:/data1 \
    --volume /data2:/data2 \
    --volume /data03:/data03 \
    --volume /data04:/data04 \
    --volume "${HOME}/.ssh:/home/docker/.ssh" \
    "${IMG}" bash -c '
        git pull --quiet && \
        pip install \
            --quiet \
            --requirement "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/jupyter-deploy/requirements.txt" && \
        pip install \
            --quiet \
            --requirement "https://raw.githubusercontent.com/ivasilyev/biopipelines-docker/master/jupyter-deploy/requirements-linux.txt"; \
        echo Web access: \"http://${IP_ADDRESS}:${PORT}/?token=${TOKEN}\";
        jupyter lab \
            --ip=0.0.0.0 \
            --no-browser \
            --NotebookApp.token=${TOKEN} \
            --port=${PORT}
    '

# Web access: `http://<ip>:31522/?token=TOKEN`
