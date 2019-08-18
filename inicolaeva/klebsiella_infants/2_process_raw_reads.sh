#!/usr/bin/env bash

# Generate sampledata
head -n 2 /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/raw_reads_pipeline.sampledata > \
    /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/raw_reads_pipeline.1.sampledata

# Clean output directory
rm -rf /data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline/*

# Copy the script file

# Run pipeline for a single sample
python3 \
    /data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline/pipeline_handler.py \
        -i /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/raw_reads_pipeline.1.sampledata \
        -o /data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline

# Run pipeline for all samples
python3 \
    /data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline/pipeline_handler.py \
        -i /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/raw_reads_pipeline.sampledata \
        -o /data1/bio/projects/inicolaeva/klebsiella_infants/pipeline

# Combine MLST results
mlst=$(find /data1/bio/projects/inicolaeva/klebsiella_infants/pipeline/11_srst2 -name *__results.txt -print)
cat <(echo "${mlst}" | head -n 1 | xargs head -qn 1) <(echo "${mlst}" | xargs tail -qn 1 | sort) | \
    tee /data1/bio/projects/inicolaeva/klebsiella_infants/pipeline/mlst.tsv

# Get help string from the main pipeline
python3 /data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline/pipeline_handler.py -h

# Run pipeline for all samples only for OrthoMCL step
python3 \
    /data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline/pipeline_handler.py \
        -i /data1/bio/projects/inicolaeva/klebsiella_infants/sample_data/raw_reads_pipeline.sampledata \
        -o /data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline \
        -s 11

# Test QUAST output
export IMG=quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0 && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} \
bash -c \
    'quast \
        -o /data1/bio/projects/inicolaeva/klebsiella_infants/test/quast \
        -t 8 \
        /data1/bio/projects/inicolaeva/klebsiella_infants/test/pipeline/04_spades/Kleb102/genome/contigs.fasta
     chmod -R 777 /data1/bio/projects/inicolaeva/klebsiella_infants/test/quast'
