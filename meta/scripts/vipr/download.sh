#!/usr/bin/env bash

# Unfortunately, the version has to be edited manually because of the reference fetch complexity.
# It may be obtained from the very last line of 'blastReadMeDownload.spg', the link is below.
# Some cookie hacks are required to obtain the file in automated order.
# Furthermore, the compressed sequence files are huge.
# Lots of time should be spent to download them via the regular fetchers, like `curl` or `wget`.
# However, a multithreaded file download & postprocessing can be safely applied it this case.

VERSION="2021.16.11"

export OUT_DIR="/data/reference/ViPR/vipr_v${VERSION}"

export IMG="p3terx/aria2-pro:latest" && \
docker pull ${IMG} && \
docker run \
    --env OUT_DIR="${OUT_DIR}" \
    --interactive \
    --net=host \
    --rm \
    --tty \
    --volume "/data:/data" \
    --volume "/data1:/data1" \
    --volume "/data2:/data2" \
    ${IMG} \
        bash -c '
            mkdir -p "${OUT_DIR}";
            for DL_URL in \
                "https://www.viprbrc.org/brcDocs/datafiles/blast/DB_new_format/NONFLU_All.nt.tar.gz" \
                "https://www.viprbrc.org/brcDocs/datafiles/blast/DB_new_format/NONFLU_All_cds.nt.tar.gz" \
                "https://www.viprbrc.org/brcDocs/datafiles/blast/DB_new_format/NONFLU_All.aa.tar.gz" \
                "https://www.viprbrc.org/brc/blastReadMeDownload.spg";
                do
                    aria2c \
                        --continue \
                        --dir="${OUT_DIR}" \
                        --file-allocation=none \
                        --max-connection-per-server=$(nproc) \
                        --split=$(nproc) \
                        "${DL_URL}";
                done;
            chmod -R a+rw "${OUT_DIR}";
        '

find "${OUT_DIR}" -name '*.tar.gz' -type f | xargs -P $(nproc) -I "{}" \
    bash -c '
        EXTRACT_DIR="${OUT_DIR}/$(echo "$(basename "{}")" | cut -d. -f1-2)";
        echo "Decompress {} into ${EXTRACT_DIR}";
        mkdir -p "${EXTRACT_DIR}";
        tar -C "${EXTRACT_DIR}" -xf "{}";
        rm -f "{}";
        echo "Processed {}";
    '
