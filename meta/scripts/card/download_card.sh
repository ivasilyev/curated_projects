#!/usr/bin/env bash

export REFERENCE_VERSION="$(
    curl -s "https://card.mcmaster.ca/download" \
    | grep 'td class=\"hidden-xs\"' \
    | grep -Eo '[0-9]+\.[0-9]+\.[0-9]+' \
    | head -n 1
)"
echo "Found CARD version: ${REFERENCE_VERSION}"

export REFERENCE_DIR="/data/reference/CARD/card_v${REFERENCE_VERSION}/"
export REFERENCE_NFASTA="${REFERENCE_DIR}card_v${REFERENCE_VERSION}.fasta"

mkdir -p "${REFERENCE_DIR}"
cd "${REFERENCE_DIR}" || exit 1

for i in \
    "https://card.mcmaster.ca/latest/data" \
    "https://card.mcmaster.ca/latest/ontology" \
    "https://card.mcmaster.ca/latest/variants"
    do
    echo "Download $(basename "${i}")"
    curl -fsSL "${i}" | tar jxf -
    done

zcat "${REFERENCE_DIR}nucleotide_fasta_protein_homolog_model_variants.fasta.gz" \
    > "${REFERENCE_NFASTA}"

export IMG="ivasilyev/bwt_filtering_pipeline_worker:latest" && \
docker pull "${IMG}" && \
docker run \
    -e REFERENCE_DIR="${REFERENCE_DIR}" \
    -e REFERENCE_NFASTA="${REFERENCE_NFASTA}" \
    -it \
    --rm \
    -v "/data:/data" \
    "${IMG}" \
    bash -c '
        cd "${REFERENCE_DIR}";
        python3 "${HOME}/scripts/cook_the_reference.py" \
        --input "${REFERENCE_NFASTA}" \
        --output "${REFERENCE_DIR}index";
    '
