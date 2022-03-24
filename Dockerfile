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
# mkdir -p ".docker/pushing_images" && cd ".docker/pushing_images" && curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/Dockerfile"
# export REPO="curated_projects" && export TAG="ivasilyev/${REPO}:latest" && docker build --network=host --tag "${REPO}" . && docker tag "${REPO}" "${TAG}" && docker push "${TAG}"  && docker push "${TAG}" && rm -f Dockerfile
