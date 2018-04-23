#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Copy 'sam2coverage.py' script into '/data1/bio/projects/ndanilova/colitis_crohn/test/' dir

# CLI preparing:

export DOCKER_IMAGE_NAME=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${DOCKER_IMAGE_NAME} python3
"""

import subprocess
import multiprocessing


def external_route(*args):
    process = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    if output:
        print(output.decode("utf-8"))


def multi_core_queue(function_name, queue):
    pool = multiprocessing.Pool()
    pool.map(function_name, queue)
    pool.close()
    pool.join()


"""
# Manual checkout

export DOCKER_ALIAS='docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ivasilyev/bwt_filtering_pipeline_worker ' 

rm -f /data2/bio/Metagenomes/custom/25_ecoli_genes/Mapped_reads/136VZK_no_hg19_25_ecoli_genes.bam \
/data2/bio/Metagenomes/custom/25_ecoli_genes/Mapped_reads/136VZK_no_hg19_25_ecoli_genes_sorted.bam \
/data1/bio/projects/ndanilova/colitis_crohn/test/136VZK_no_hg19_25_ecoli_genes.tsv

${DOCKER_ALIAS} samtools view -t /data/reference/custom/25_ecoli_genes/index/25_ecoli_genes_samtools.fai \
-bS /data2/bio/Metagenomes/custom/25_ecoli_genes/Mapped_reads/136VZK_no_hg19_25_ecoli_genes.sam \
-o /data2/bio/Metagenomes/custom/25_ecoli_genes/Mapped_reads/136VZK_no_hg19_25_ecoli_genes.bam \
-@ 32

file /data2/bio/Metagenomes/custom/25_ecoli_genes/Mapped_reads/136VZK_no_hg19_25_ecoli_genes.bam

${DOCKER_ALIAS} samtools sort /data2/bio/Metagenomes/custom/25_ecoli_genes/Mapped_reads/136VZK_no_hg19_25_ecoli_genes.bam \
-o /data2/bio/Metagenomes/custom/25_ecoli_genes/Mapped_reads/136VZK_no_hg19_25_ecoli_genes_sorted.bam \
-@ 32 

file /data2/bio/Metagenomes/custom/25_ecoli_genes/Mapped_reads/136VZK_no_hg19_25_ecoli_genes_sorted.bam

${DOCKER_ALIAS} genomeCoverageBed -ibam /data2/bio/Metagenomes/custom/25_ecoli_genes/Mapped_reads/136VZK_no_hg19_25_ecoli_genes_sorted.bam \
> /data1/bio/projects/ndanilova/colitis_crohn/test/136VZK_no_hg19_25_ecoli_genes.tsv

file /data1/bio/projects/ndanilova/colitis_crohn/test/136VZK_no_hg19_25_ecoli_genes.tsv

"""


def mp_get_coverage(external_input):
    external_route('python3',  '/data1/bio/projects/ndanilova/colitis_crohn/test/sam2coverage.py',
                   '-i', external_input,
                   '-f', "/data/reference/custom/25_ecoli_genes/index/25_ecoli_genes_samtools.fai",
                   '-g', "/data/reference/custom/25_ecoli_genes/index/25_ecoli_genes_samtools.genome",
                   "-a", "/data/reference/custom/25_ecoli_genes/index/25_ecoli_genes_annotation.txt",
                   '-o', "/data2/bio/Metagenomes/custom/25_ecoli_genes")


def sam2coverage():
    [mp_get_coverage(i) for i in samFileNamesList]
    # multi_core_queue(mp_get_coverage, samFileNamesList)


samFileNamesList = sorted(subprocess.getoutput("ls -d /data2/bio/Metagenomes/custom/25_ecoli_genes/Mapped_reads/*.sam").split("\n"))
sam2coverage()

# subprocess.getoutput("samtools view -Su /data2/bio/Metagenomes/custom/25_ecoli_genes/Mapped_reads/136VZK_no_hg19_25_ecoli_genes.sam -@ 32 | \
#                      samtools sort - -o /data1/bio/projects/ndanilova/colitis_crohn/test/136VZK_no_hg19_25_ecoli_genes.sorted.bam -@ 32 ")
