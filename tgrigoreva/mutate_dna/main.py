#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -e DISPLAY=$DISPLAY -it ${DOCKER_IMAGE_NAME} python3
"""

from tgrigoreva.mutate_dna.FASTQArray import FASTAArray
from meta.scripts.Utilities import Utilities
import subprocess
import re
import gc

project_dir = ""
r1_fastq_archives_list, r2_fastq_archives_list = [sorted([i.strip() for i in subprocess.getoutput("find {} -name *R{}*.fastq.gz".format(project_dir, j)).split("\n") if len(i.strip()) > 0]) for j in (1, 2)]
max_id = max([int(re.findall("[0-9]{6}", i)[0]) for i in r1_fastq_archives_list + r2_fastq_archives_list])
output_files_list = []
for r1_fastq_archive, r2_fastq_archive in zip(r1_fastq_archives_list, r2_fastq_archives_list):
    for copy_number in [1, 2]:
        max_id += 1
        for fastq_archive in [r1_fastq_archive, r2_fastq_archive]:
            output_file = "{}{}".format(Utilities.ends_with_slash("/".join(fastq_archive.split("/")[:-1])), re.sub("\.gz$", "", re.sub("[0-9]{6}", str(max_id).zfill(6), fastq_archive.split("/")[-1])))
            print("Loading file '{}'".format(fastq_archive))
            fq_array = FASTAArray(subprocess.getoutput("zcat {}".format(fastq_archive)))
            print("Loaded file '{}'".format(fastq_archive))
            fq_array.parse_fastq(output_file)
            del fq_array
            output_files_list.append(output_file)
            print("Saved file '{}'".format(output_file))
            gc.collect()


def mp_gzip_file(file):
    print("Compressing file '{}'".format(file))
    print(subprocess.getoutput("gzip -9 -c {a} > {a}.gz".format(a=file)))
    print("Compressed file '{}'".format(file))


print("Compressing {} files".format(len(output_files_list)))
print(Utilities.multi_core_queue(mp_gzip_file, output_files_list))
