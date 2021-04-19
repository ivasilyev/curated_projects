#!/usr/bin/env bash

# The SAMPLEDATA_MASK and METADATA_TSV variables are defined externally
mkdir -p "${SAMPLEDATA_MASK}/logs"
cd "${SAMPLEDATA_MASK}" || exit 1

NPROC="$(grep -c '^processor' /proc/cpuinfo)"

# From https://antonioggsousa.github.io/tutorial/example/
echo Import and convert paired-end FASTQ files to a QIIME2 artifact
mkdir -p "demultiplexed_reads"
qiime tools import --input-format PairedEndFastqManifestPhred33 \
  --input-path "../../sample_data/chopped/qiime2_sample_data_${SAMPLEDATA_MASK}.csv" \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --output-path "demultiplexed_reads/demultiplexed_PE_reads.qza" \
  |& tee "logs/tools import.log"

echo Demultiplex and summarize sequences
mkdir -p "visualized"
qiime demux summarize --verbose \
  --i-data "demultiplexed_reads/demultiplexed_PE_reads.qza" \
  --o-visualization "visualized/demux_PE_reads.qzv" \
  |& tee "logs/demux summarize demux_PE_reads.log"

#echo Visualize the plots
#qiime tools view "visualized/demux_PE_reads.qzv"
#
#echo Quality filtering
#qiime quality-filter q-score --verbose \
#  --i-demux "demultiplexed_reads/demultiplexed_PE_reads.qza" \
#  --o-filtered-sequences "demux-paired-end-filtered.qza" \
#  --o-filter-stats "demux-filter-stats.qza"

#echo Export back the reads for a detailed check
#qiime tools export \
#  --input-path "demultiplexed_reads/demultiplexed_PE_reads.qza" \
#  --output-path "import-check"

echo DADA2 denoising
mkdir -p "dada2"
qiime dada2 denoise-paired --verbose \
  --p-trunc-len-f 225 --p-trunc-len-r 225 --p-n-reads-learn 30000 \
  --p-n-threads "${NPROC}" \
  --i-demultiplexed-seqs "demultiplexed_reads/demultiplexed_PE_reads.qza" \
  --o-representative-sequences "dada2/dada2_representative_sequences.qza" \
  --o-table "dada2/dada2_frequency_table.qza" \
  --o-denoising-stats "dada2/dada2_denoising_statistics.qza" \
  |& tee "logs/dada2 denoise-paired.log"

echo Summarize stats
qiime feature-table summarize --verbose \
  --i-table "dada2/dada2_frequency_table.qza"\
  --o-visualization "visualized/table-dada2.qzv"
qiime feature-table tabulate-seqs --verbose \
  --i-data "dada2/dada2_representative_sequences.qza" \
  --o-visualization "visualized/rep-seqs-dada2.qzv"

echo Assign taxonomy
mkdir -p "classified"
# --p-n-jobs, The maximum number of concurrently worker processes. If -1 all CPUs are used. If 1 is given, no parallel computing code is used at all, which is useful for debugging. For n-jobs below -1, (n_cpus + 1 + n-jobs) are used. Thus for n-jobs = -2, all CPUs but one are used.
qiime feature-classifier classify-sklearn --verbose \
  --p-n-jobs "-1" \
  --p-reads-per-batch 10000 \
  --i-classifier "/data/reference/SILVA/SILVA_v138/SILVA-138-SSURef-full-length-classifier.qza" \
  --i-reads "dada2/dada2_representative_sequences.qza" \
  --o-classification "classified/classified_taxonomy.qza" \
  |& tee "logs/feature-classifier classify-sklearn.log"

echo Make an Amplicon Sequence Variant table
qiime metadata tabulate --verbose \
  --m-input-file "classified/classified_taxonomy.qza" \
  --o-visualization "visualized/classified_taxonomy.qzv" \
  |& tee "logs/metadata tabulate classified_taxonomy.log"

echo Make a prokaryotic profile
qiime taxa barplot --verbose \
  --m-metadata-file "${METADATA_TSV}" \
  --i-table "dada2/dada2_frequency_table.qza"\
  --i-taxonomy "classified/classified_taxonomy.qza" \
  --o-visualization "visualized/taxonomy_barplots.qzv" \
  |& tee "logs/taxa barplot.log"

echo Perform de novo multiple sequence alignment
mkdir -p "aligned"
qiime alignment mafft --verbose \
  --p-n-threads "${NPROC}" \
  --i-sequences "dada2/dada2_representative_sequences.qza" \
  --o-alignment "aligned/aligned_sequences.qza" \
  |& tee "logs/alignment mafft.log"

echo Filter the unconserved and highly variable and gapped columns to avoid overestimate distances
mkdir -p "masked"
qiime alignment mask --verbose \
  --i-alignment "aligned/aligned_sequences.qza" \
  --o-masked-alignment "masked/masked_aligned_sequences.qza" \
  |& tee "logs/alignment mask.log"

echo Build a phylogenetic ML tree
mkdir -p "unrooted_trees"
qiime phylogeny fasttree --verbose \
  --p-n-threads "${NPROC}" \
  --i-alignment "masked/masked_aligned_sequences.qza" \
  --o-tree "unrooted_trees/unrooted_tree.qza" \
  |& tee "logs/phylogeny fasttree.log"

echo Root the unrooted tree based on the midpoint rooting method
mkdir -p "rooted_trees"
qiime phylogeny midpoint-root --verbose \
  --i-tree "unrooted_trees/unrooted_tree.qza" \
  --o-rooted-tree "rooted_trees/rooted_tree.qza" \
  |& tee "logs/phylogeny midpoint-root.log"

echo Analyze the core diversity using the phylogenetic pipeline
qiime diversity core-metrics-phylogenetic --verbose --p-sampling-depth 1000 \
  --p-n-jobs-or-threads "${NPROC}" \
  --m-metadata-file "${METADATA_TSV}" \
  --i-phylogeny "rooted_trees/rooted_tree.qza" \
  --i-table "dada2/dada2_frequency_table.qza"\
  --output-dir "phylogenetic_core_metrics" \
  |& tee "logs/diversity core-metrics-phylogenetic.log"

echo Build rarefaction plots
qiime diversity alpha-rarefaction --verbose --p-max-depth 2000 \
  --m-metadata-file "${METADATA_TSV}" \
  --i-table "dada2/dada2_frequency_table.qza" \
  --i-phylogeny "rooted_trees/rooted_tree.qza" \
  --o-visualization "visualized/alpha_rarefaction.qzv" \
  |& tee "logs/diversity alpha-rarefaction.log"

echo Join paired-end reads
# Does not scale much past 4 threads.
mkdir -p "joined_reads"
qiime vsearch join-pairs --verbose --p-allowmergestagger \
  --p-threads 8 \
  --i-demultiplexed-seqs "demultiplexed_reads/demultiplexed_PE_reads.qza" \
  --o-joined-sequences "joined_reads/joined_PE_reads.qza" \
  |& tee "logs/vsearch join-pairs.log"

echo Filter based on Q scores
mkdir -p "q_score_filter"
qiime quality-filter q-score --verbose \
  --i-demux "joined_reads/joined_PE_reads.qza" \
  --o-filtered-sequences "q_score_filter/sequences_filtered_by_q_score.qza" \
  --o-filter-stats "q_score_filter/filtering_statistics.qza" \
  |& tee "logs/quality-filter q-score.log"

echo Dereplicate sequences
mkdir -p "dereplicated"
qiime vsearch dereplicate-sequences --verbose \
  --i-sequences "q_score_filter/sequences_filtered_by_q_score.qza" \
  --o-dereplicated-table "dereplicated/dereplicated_frequency_table.qza" \
  --o-dereplicated-sequences "dereplicated/dereplicated_sequences.qza" \
  |& tee "logs/vsearch dereplicate-sequences.log"

echo Cluster closed references at 97%
mkdir -p "closed_reference"
qiime vsearch cluster-features-closed-reference --verbose --p-perc-identity 0.97 \
  --p-threads "${NPROC}" \
  --i-reference-sequences "/data/reference/SILVA/SILVA_v138/SILVA-138-SSURef-Full-Seqs.qza" \
  --i-table "dereplicated/dereplicated_frequency_table.qza" \
  --i-sequences "dereplicated/dereplicated_sequences.qza" \
  --o-clustered-table "closed_reference/closed_reference_clustered_table.qza" \
  --o-clustered-sequences "closed_reference/closed_reference_clustered_sequences.qza" \
  --o-unmatched-sequences "closed_reference/closed_reference_unmatched_sequences.qza" \
  |& tee "logs/vsearch cluster-features-closed-reference.log"

echo Export an OTU table
mkdir -p "biom"
qiime tools export \
  --input-path "closed_reference/closed_reference_clustered_table.qza" \
  --output-path "biom" \
  |& tee "logs/tools export base biom.log"
# Output: 'feature-table.biom'

echo Convert biom to JSON
biom convert --to-json \
  --input-fp "biom/feature-table.biom" \
  --output-fp "biom/OTUs.json" \
  |& tee "logs/biom convert base json.log"

echo Convert biom to TSV
biom convert --to-tsv \
  --input-fp "biom/feature-table.biom" \
  --output-fp "biom/OTUs.tsv" \
  |& tee "logs/biom convert base tsv.log"

echo Create table with taxonomy annotations
qiime tools export \
  --input-path "classified/classified_taxonomy.qza" \
  --output-path "biom" \
  |& tee "logs/tools export base tsv.log"
# Output: 'taxonomy.tsv'

echo Annotate biom with taxonomy data
< "biom/taxonomy.tsv" \
  sed 's|Feature ID\tTaxon\tConfidence|#OTUID\ttaxonomy\tconfidence|' \
  > "biom/taxa.tsv"

biom add-metadata \
  --sc-separated "taxonomy" \
  --observation-metadata-fp "biom/taxa.tsv" \
  --input-fp "biom/feature-table.biom" \
  --output-fp "biom/OTUs_with_taxa.biom" \
  |& tee "logs/biom add-metadata.log"

echo Export the aligned sequences
qiime tools export \
  --input-path "closed_reference/closed_reference_clustered_sequences.qza" \
  --output-path "closed_reference" \
  |& tee "logs/tools export fasta.log"
# Output: 'dna-sequences.fasta'

chmod -R 777 "$(pwd)"
cd ..
exit 0
