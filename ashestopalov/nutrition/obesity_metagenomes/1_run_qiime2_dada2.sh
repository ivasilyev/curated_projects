#!/usr/bin/env bash

mkdir -p "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/qiime2"

export IMG=qiime2/core:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 --net=host -it ${IMG} bash

cd "/data1/bio/projects/ashestopalov/nutrition/obesity_metagenomes/qiime2"


# node5
export SSRC="stool"
# node6
export SSRC="blood"

mkdir -p ${SSRC}
chmod -R 777 ${SSRC}
cd ${SSRC}

export SMETADATA="../../sample_data/qiime2_meta_data_${SSRC}.tsv"

# From https://antonioggsousa.github.io/tutorial/example/
echo Import and convert fastq files to QIIME2 artifact
qiime tools import --input-format PairedEndFastqManifestPhred33 \
  --input-path "../../sample_data/qiime2_sample_data_${SSRC}.csv" \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --output-path "demux-paired-end.qza"

echo Demultiplex and summarize sequences
qiime demux summarize --verbose \
  --i-data "demux-paired-end.qza" \
  --o-visualization "demux-paired-end.qzv"

#echo Visualize the plots
#qiime tools view "demux-paired-end.qzv"
#
#echo Quality filtering
#qiime quality-filter q-score --verbose \
#  --i-demux "demux-paired-end.qza" \
#  --o-filtered-sequences "demux-paired-end-filtered.qza" \
#  --o-filter-stats "demux-filter-stats.qza"

#echo Export back the reads for a detailed check
#qiime tools export \
#  --input-path "demux-paired-end.qza" \
#  --output-path "import-check"

echo DADA2 denoising
qiime dada2 denoise-paired --verbose \
  --p-trunc-len-f 225 --p-trunc-len-r 225 --p-n-reads-learn 30000 --p-n-threads 0 \
  --i-demultiplexed-seqs "demux-paired-end.qza" \
  --o-representative-sequences "rep-seqs-dada2.qza" \
  --o-table "table-dada2.qza" \
  --output-dir "dada2-dmx-pe"

#echo Summarize stats
#qiime feature-table summarize --verbose \
#  --i-table "table-dada2.qza" \
#  --o-visualization "table-dada2.qzv"
#qiime feature-table tabulate-seqs --verbose \
#  --i-data "rep-seqs-dada2.qza" \
#  --o-visualization "rep-seqs-dada2.qzv"

echo Assign taxonomy
qiime feature-classifier classify-sklearn --verbose \
  --i-classifier "/data/reference/SILVA/SILVA_v138/SILVA-138-SSURef-full-length-classifier.qza" \
  --i-reads "rep-seqs-dada2.qza" \
  --o-classification "taxonomy-rep-seqs-dada2.qza"

echo Make an Amplicon Sequence Variant table
qiime metadata tabulate \
  --m-input-file "taxonomy-rep-seqs-dada2.qza" \
  --o-visualization "taxonomy-rep-seqs-dada2.qzv "

echo Make a prokaryotic profile
qiime taxa barplot \
  --m-metadata-file ${SMETADATA} \
  --i-table "table-dada2.qza" \
  --i-taxonomy "taxonomy-rep-seqs-dada2.qza" \
  --o-visualization "taxonomy-bar-plots.qzv"

echo Align multi sequences
qiime alignment mafft \
  --i-sequences "rep-seqs-dada2.qza" \
  --o-alignment "mafft-rep-seqs-dada2.qza"

echo Eliminate the highly variable positions to avoid overestimate distances
qiime alignment mask \
  --i-alignment "mafft-rep-seqs-dada2.qza" \
  --o-masked-alignment "masked-msa-rep-seqs-dada2.qza"

echo Build a ML tree
qiime phylogeny fasttree \
  --i-alignment "masked-msa-rep-seqs-dada2.qza" \
  --o-tree "unroot-ml-tree-masked.qza"

echo Root the unrooted tree based on the midpoint rooting method
qiime phylogeny midpoint-root \
  --i-tree "unroot-ml-tree-masked.qza" \
  --o-rooted-tree "root-ml-tree.qza"

echo Analyze the core diversity using the phylogenetic pipeline
qiime diversity core-metrics-phylogenetic --p-sampling-depth 20000 \
  --m-metadata-file ${SMETADATA} \
  --i-phylogeny "root-ml-tree.qza" \
  --i-table "table-dada2.qza" \
  --output-dir "core-metrics-results"

echo Build rarefaction plots
qiime diversity alpha-rarefaction --p-max-depth 20000 \
  --m-metadata-file ${SMETADATA} \
  --i-table "table-dada2.qza" \
  --i-phylogeny "root-ml-tree.qza" \
  --o-visualization "alpha-rarefaction.qzv"

echo Join paired-end reads
qiime vsearch join-pairs --p-allowmergestagger \
  --i-demultiplexed-seqs "demux-paired-end.qza" \
  --o-joined-sequences "dmx-jpe.qza"

echo Filter based on Q scores
qiime quality-filter q-score \
  --i-demux "dmx-jpe.qza" \
  --o-filtered-sequences "dmx-jpe-filter.qza" \
  --o-filter-stats "dmx-jpe-filter-stats.qza"

echo Dereplicate sequences
qiime vsearch dereplicate-sequences \
  --i-sequences "dmx-jpe-filter.qza" \
  --o-dereplicated-table "drpl-tbl.qza" \
  --o-dereplicated-sequences "drpl-seqs.qza"

echo Cluster closed references at 97%
qiime vsearch cluster-features-closed-reference --p-perc-identity 0.97 \
  --i-reference-sequences "/data/reference/SILVA/SILVA_v138/SILVA-138-SSURef-Full-Seqs.qza" \
  --i-table "drpl-tbl.qza" \
  --i-sequences "drpl-seqs.qza" \
  --o-clustered-table "tbl-cr-97.qza" \
  --o-clustered-sequences "rep-seqs-cr-97.qza" \
  --o-unmatched-sequences "unmatched-cr-97.qza"

echo Export the OTU table
qiime tools export \
  --input-path "tbl-cr-97.qza" \
  --output-path .

echo Convert biom to json
biom convert --to-json \
  --input-fp "feature-table.biom" \
  --output-fp "feature-table.json"

echo Convert biom to TSV
biom convert --to-tsv \
  --input-fp "feature-table.biom" \
  --output-fp "meta-predic.tsv"

echo Export the aligned sequences
mkdir -p "cr-97"
qiime tools export \
  --input-path "rep-seqs-cr-97.qza" \
  --output-path "cr-97"
