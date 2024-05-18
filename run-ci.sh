#!/usr/bin/env bash

export IMG_NAME="ivasilyev/curated_projects"
export IMG_TAG="latest"
export SLEEP_TIME="10h"

export IMG_DIR="/tmp/${IMG_NAME}/"
export COMMIT_FILE="${IMG_DIR}commit"

mkdir -p "${IMG_DIR}"

cd "${IMG_DIR}"

export COMMIT="$(cat "${COMMIT_FILE}")"

export LATEST_COMMIT="$(
    curl -fsSL \
        "https://api.github.com/repos/ivasilyev/curated_projects/commits?path=meta&page=1&per_page=1" \
    | grep "comments_url" \
    | grep \
        --only-matching \
        --perl-regexp \
        '(?<=commits/)[^/]*(?=/comments)'
)"

if [ "${COMMIT}" == "${LATEST_COMMIT}" ]
    then
    echo "Already built commit: ${LATEST_COMMIT}"
else
    echo "Build ${COMMIT} -> ${LATEST_COMMIT}"

    curl -fsSLO "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/Dockerfile"

    docker build \
        --network host \
        --tag "${IMG_NAME}" \
        "${IMG_DIR}" && \
    docker tag \
        "${IMG_NAME}" \
        "${IMG_TAG}" && \
    docker push "${IMG_NAME}:${IMG_TAG}"

    rm -rf ./*

    echo "${LATEST_COMMIT}" > "${COMMIT_FILE}"

fi

echo "Pause after ${LATEST_COMMIT}"
sleep "${SLEEP_TIME}"
