FROM ivasilyev/jupyter-deploy:latest

RUN mkdir -p /home/docker/scripts && \
    cd /home/docker/scripts/ && \
    git clone https://github.com/ivasilyev/curated_projects.git && \
    git clone https://github.com/ivasilyev/statistical_tools.git

WORKDIR /home/docker/scripts/curated_projects/

ENV PYTHONPATH=/home/docker/scripts/curated_projects/
ENV UTILS_DIR=/home/docker/scripts/curated_projects/meta/scripts/

CMD ["/bin/bash"]

# MANUAL BUILD COMMAND:
# export REPO="curated_projects" && export TAG="ivasilyev/${REPO}:latest" && docker build --network=host --tag "${REPO}" . && docker tag "${REPO}" "${TAG}" && docker push "${TAG}"
