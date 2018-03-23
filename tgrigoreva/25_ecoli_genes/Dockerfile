FROM debian:jessie

# CLI: docker pull debian:jessie && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it debian:jessie bash

RUN apt-get -y update && \
    apt-get -y install \
            ca-certificates \
            curl \
            python3 \
            python3-pip && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip3 install bs4 \
                 requests

CMD ["/bin/bash"]

# MANUAL BUILD COMMAND:
# docker build -t env_25_ecoli_genes . && docker tag env_25_ecoli_genes ivasilyev/env_25_ecoli_genes:latest && docker push ivasilyev/env_25_ecoli_genes:latest
