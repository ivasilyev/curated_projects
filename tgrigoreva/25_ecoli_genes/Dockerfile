FROM debian:jessie

# CLI: docker pull debian:jessie && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it debian:jessie bash

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update && \
    apt-get -y install \
            # Base packets
            ca-certificates \
            curl \
            python3 \
            python3-pip \
            # Required by lxml
            libxml2-dev \
            libxslt1-dev \
            zlib1g-dev \
    && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip3 install bs4 \
                 lxml \
                 requests

CMD ["/bin/bash"]

# MANUAL BUILD COMMAND:
# docker build -t env_25_ecoli_genes . && docker tag env_25_ecoli_genes ivasilyev/env_25_ecoli_genes:latest && docker push ivasilyev/env_25_ecoli_genes:latest
