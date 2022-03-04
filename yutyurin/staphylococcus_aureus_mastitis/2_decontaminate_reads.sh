#!/usr/bin/env bash

DECONTAMINATION_DIR="${ROOT_DIR}decontamination/"
mkdir -p "${DECONTAMINATION_DIR}"
cd "${DECONTAMINATION_DIR}" || exit 1

for i in ${!SAMPLE_NAMES[@]}
    do
    SAMPLE_NAME="${SAMPLE_NAMES[$i]}"

    export IMG="ivasilyev/curated_projects:latest" && \
    docker pull "${IMG}" && \
    docker run \
        --env IN_FNA="${ROOT_DIR}pga-pe-pipeline/06_plasmid_merger/${SAMPLE_NAME}/${SAMPLE_NAME}_genome.fna" \
        --env CONTAMINATION_REPORT="${DECONTAMINATION_DIR}Contamination_${SAMPLE_NAME}_genome.txt" \
        --env OUT_FNA="${DECONTAMINATION_DIR}${SAMPLE_NAME}_genome.fna" \
        --interactive \
        --net=host \
        --rm \
        --tty \
        --volume /data:/data \
        --volume /data1:/data1 \
        --volume /data2:/data2 \
        --volume /data03:/data03 \
        --volume /data04:/data04 \
        ${IMG} \
            bash -c '
                git pull --quiet && \
                python3 "${HOME}/scripts/curated_projects/meta/scripts/ncbi_contamination_remover.py" \
                    --input "${IN_FNA}" \
                    --contamination "${CONTAMINATION_REPORT}" \
                    --output "${OUT_FNA}"
            '
    done
