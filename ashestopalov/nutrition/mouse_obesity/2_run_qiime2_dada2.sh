#!/usr/bin/env bash

function log {
    echo "[$(date '+%d-%m-%Y %H:%M:%S')] $@"
}

# QIIME2_DIR, SAMPLEDATA_CSV and METADATA_TSV variables are defined externally
log "Run QIIME2 in ${QIIME2_DIR}"

export LOG_DIR="${QIIME2_DIR}logs/"
export CONSENSUS_THRESHOLD=97
export NPROC="$(grep -c '^processor' "/proc/cpuinfo")"
export TAXA_REFERENCE_FEATURES="/data/reference/SILVA/SILVA_v138/Silva-v138-full-length-seq-taxonomy.qza"
export TAXA_REFERENCE_CLASSIFIER="/data/reference/SILVA/SILVA_v138/SILVA-138-SSURef-full-length-classifier.qza"
export TAXA_REFERENCE_SEQUENCES="/data/reference/SILVA/SILVA_v138/SILVA-138-SSURef-Full-Seqs.qza"
export TAXA_REFERENCE_HEADER="/data/reference/SILVA/SILVA_v138/SILVA_138_Taxonomy_headed.tsv"


mkdir -p "${LOG_DIR}"
cd "${QIIME2_DIR}" || exit 1

log Import and convert pre-demultiplexed paired-end FASTQ files to QIIME2 artifact
mkdir -p "demultiplexed_reads/"
qiime tools import --input-format PairedEndFastqManifestPhred33 \
    --input-path "${SAMPLEDATA_CSV}" \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --output-path "demultiplexed_reads/demultiplexed_PE_reads.qza" \
    |& tee "${LOG_DIR}tools import.log"

log Summarize sequences
mkdir -p "visualizations/"
qiime demux summarize \
    --i-data "demultiplexed_reads/demultiplexed_PE_reads.qza" \
    --o-visualization "visualizations/demultiplexed_PE_reads.qzv" \
     --verbose \
    |& tee "${LOG_DIR}demux summarize demux_PE_reads.log"

log DADA2 denoising
mkdir -p "dada2/"
qiime dada2 denoise-paired \
    --p-trunc-len-f 225 \
    --p-trunc-len-r 225 \
    --p-n-reads-learn 30000 \
    --p-n-threads "${NPROC}" \
    --i-demultiplexed-seqs "demultiplexed_reads/demultiplexed_PE_reads.qza" \
    --o-representative-sequences "dada2/dada2_representative_sequences.qza" \
    --o-table "dada2/dada2_frequency_table.qza" \
    --o-denoising-stats "dada2/dada2_denoising_statistics.qza" \
    --verbose \
    |& tee "${LOG_DIR}dada2 denoise-paired.log"
qiime metadata tabulate \
    --m-input-file "dada2/dada2_denoising_statistics.qza" \
    --o-visualization "visualizations/dada2_denoising_statistics.qza" \
    --verbose \
    |& tee "${LOG_DIR} metadata tabulate.log"

log Summarize statistics
qiime feature_table summarize \
    --i-table "dada2/dada2_frequency_table.qza"\
    --o-visualization "visualizations/dada2_frequency_table.qzv" \
    --m-sample-metadata-file "${METADATA_TSV}" \
    --verbose
qiime feature_table tabulate-seqs \
    --i-data "dada2/dada2_representative_sequences.qza" \
    --o-visualization "visualizations/dada2_representative_sequences.qzv" \
    --verbose

log Assign taxonomy
mkdir -p "taxonomy/"
# --p-n-jobs, The maximum number of concurrently worker processes. If -1 all CPUs are used. If 1 is given, no parallel computing code is used at all, which is useful for debugging. For n-jobs below -1, (n_cpus + 1 + n-jobs) are used. Thus for n-jobs = -2, all CPUs but one are used.
qiime feature-classifier classify-sklearn \
    --p-n-jobs "-1" \
    --p-reads-per-batch 10000 \
    --i-classifier "${TAXA_REFERENCE_CLASSIFIER}" \
    --i-reads "dada2/dada2_representative_sequences.qza" \
    --o-classification "taxonomy/classified_taxonomy.qza" \
    --verbose \
    |& tee "${LOG_DIR}feature-classifier classify-sklearn.log"

log Create Amplicon Sequence Variant \(ASV\) table
qiime metadata tabulate \
    --m-input-file "taxonomy/classified_taxonomy.qza" \
    --o-visualization "visualizations/classified_taxonomy.qzv" \
    --verbose \
    |& tee "${LOG_DIR}metadata tabulate classified_taxonomy.log"

log Make prokaryotic profile
qiime taxa barplot \
    --m-metadata-file "${METADATA_TSV}" \
    --i-table "dada2/dada2_frequency_table.qza" \
    --i-taxonomy "taxonomy/classified_taxonomy.qza" \
    --o-visualization "visualizations/taxonomy_barplots.qzv" \
    --verbose \
    |& tee "${LOG_DIR}taxa barplot.log"

log Join paired-end reads
# Threads number must be within [0, 8].
mkdir -p "joined_reads/"
qiime vsearch join-pairs \
    --p-threads 8 \
    --i-demultiplexed-seqs "demultiplexed_reads/demultiplexed_PE_reads.qza" \
    --p-allowmergestagger \
    --o-joined-sequences "joined_reads/joined_PE_reads.qza" \
    --verbose \
    |& tee "${LOG_DIR}vsearch join-pairs.log"

log Filter based on Q scores
mkdir -p "q_score_filtered_reads/"
qiime quality-filter q-score \
    --i-demux "joined_reads/joined_PE_reads.qza" \
    --o-filtered-sequences "q_score_filtered_reads/sequences_filtered_by_q_score.qza" \
    --o-filter-stats "q_score_filtered_reads/filtering_statistics.qza" \
    --verbose \
    |& tee "${LOG_DIR}quality-filter q-score.log"

log Dereplicate sequences
mkdir -p "dereplicated/"
qiime vsearch dereplicate-sequences \
    --i-sequences "q_score_filtered_reads/sequences_filtered_by_q_score.qza" \
    --o-dereplicated-table "dereplicated/dereplicated_frequency_table.qza" \
    --o-dereplicated-sequences "dereplicated/dereplicated_sequences.qza" \
    --verbose \
    |& tee "${LOG_DIR}vsearch dereplicate-sequences.log"

log Cluster closed references at ${CONSENSUS_THRESHOLD}%
mkdir -p "closed_references/"
qiime vsearch cluster-features-closed-reference \
    --p-threads "${NPROC}" \
    --i-reference-sequences "${TAXA_REFERENCE_SEQUENCES}" \
    --i-table "dereplicated/dereplicated_frequency_table.qza" \
    --i-sequences "dereplicated/dereplicated_sequences.qza" \
    --o-clustered-table "closed_references/closed_reference_clustered_table.qza" \
    --o-clustered-sequences "closed_references/closed_reference_clustered_sequences.qza" \
    --o-unmatched-sequences "closed_references/closed_reference_unmatched_sequences.qza" \
    --p-perc-identity 0.${CONSENSUS_THRESHOLD} \
    --verbose \
    |& tee "${LOG_DIR}vsearch cluster-features-closed-reference.log"

log Export the aligned sequences
qiime tools export \
    --input-path "closed_references/closed_reference_clustered_sequences.qza" \
    --output-format DNASequencesDirectoryFormat \
    --output-path "closed_reference" \
    |& tee "${LOG_DIR}tools export fasta.log"
# Output: 'dna-sequences.fasta'

log Export an ASV table
mkdir -p "bioms/"
qiime tools export \
    --input-path "closed_references/closed_reference_clustered_table.qza" \
    --output-path "bioms/" \
    --output-format BIOMV210DirFmt \
    |& tee "${LOG_DIR}tools export feature-table.biom.log"
# Output: 'feature-table.biom'

log Annotate biom with taxonomy data
biom add-metadata \
    --sc-separated "taxonomy" \
    --observation-metadata-fp "${TAXA_REFERENCE_HEADER}" \
    --sample-metadata-fp "${METADATA_TSV}" \
    --input-fp "bioms/feature-table.biom" \
    --output-fp "bioms/ASVs_with_taxa.biom" \
    |& tee "${LOG_DIR}biom add-metadata.log"

log Convert biom to JSON
biom convert \
    --to-json \
    --input-fp "bioms/ASVs_with_taxa.biom" \
    --output-fp "bioms/ASVs_with_taxa.json" \
    |& tee "${LOG_DIR}biom convert json.log"

log Convert biom to TSV
biom convert \
    --to-tsv \
    --input-fp "bioms/ASVs_with_taxa.biom" \
    --output-fp "bioms/ASVs_with_taxa.tsv" \
    --header-key "taxonomy" \
    |& tee "${LOG_DIR}biom convert taxa tsv.log"



log Perform de novo multiple sequence alignment
mkdir -p "alignments/"
qiime alignment mafft \
    --p-n-threads "${NPROC}" \
    --i-sequences "dada2/dada2_representative_sequences.qza" \
    --o-alignment "alignments/aligned_sequences.qza" \
    --verbose \
    |& tee "${LOG_DIR}alignment mafft.log"

log Filter the unconserved and highly variable and gapped columns to avoid overestimate distances
mkdir -p "masked_alignments/"
qiime alignment mask \
    --i-alignment "alignments/aligned_sequences.qza" \
    --o-masked-alignment "masked_alignments/masked_aligned_sequences.qza" \
    --verbose \
    |& tee "${LOG_DIR}alignment mask.log"

log Build a phylogenetic ML tree
mkdir -p "unrooted_trees/"
qiime phylogeny fasttree \
    --p-n-threads "${NPROC}" \
    --i-alignment "masked_alignments/masked_aligned_sequences.qza" \
    --o-tree "unrooted_trees/unrooted_tree.qza" \
    --verbose \
    |& tee "${LOG_DIR}phylogeny fasttree.log"

log Root the unrooted tree based on the midpoint rooting method
mkdir -p "rooted_trees/"
qiime phylogeny midpoint-root \
    --i-tree "unrooted_trees/unrooted_tree.qza" \
    --o-rooted-tree "rooted_trees/rooted_tree.qza" \
    --verbose \
    |& tee "${LOG_DIR}phylogeny midpoint-root.log"

qiime tools export \
    --input-path "dada2/dada2_frequency_table.qza" \
    --output-format BIOMV210DirFmt \
    --output-path "dada2/"
# Output: 'feature-table.biom'
biom convert \
    --to-tsv \
    --input-fp "dada2/feature-table.biom" \
    --output-fp "dada2/feature-table.tsv" \
    |& tee "${LOG_DIR}biom convert dada2 tsv.log"


# The first 2 lines are '# Constructed from biom file' and header
export DENOISED_SAMPLES=$(( $(wc -l "dada2/feature-table.tsv" | awk '{ print $1 }') - 2 ))

log Analyze the core diversity using the phylogenetic pipeline
# '--output-dir' must not exist!
rm -rf "phylogenetic_core_metrics/"
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "rooted_trees/rooted_tree.qza" \
    --i-table "dada2/dada2_frequency_table.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --output-dir "phylogenetic_core_metrics/" \
    --p-n-jobs-or-threads "${NPROC}" \
    --p-sampling-depth 10 \
    --verbose \
    |& tee "${LOG_DIR}diversity core-metrics-phylogenetic.log"
qiime metadata tabulate \
    --m-input-file "phylogenetic_core_metrics/faith_pd_vector.qza" \
    --o-visualization "visualizations/faith-pd-group-significance.qzv" \
    --verbose
qiime tools export \
    --input-path "phylogenetic_core_metrics/faith_pd_vector.qza" \
    --output-format AlphaDiversityDirectoryFormat \
    --output-path "phylogenetic_core_metrics/" \
    |& tee "${LOG_DIR}tools export faith_pd_vector.log"
# Output: alpha-diversity.tsv

log Visualize alpha diversity
qiime diversity alpha-group-significance \
    --i-alpha-diversity "phylogenetic_core_metrics/faith_pd_vector.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --o-visualization "visualizations/alpha_faith_pd_group_significance.qzv" \
    --verbose \
    |& tee "${LOG_DIR}diversity alpha-group-significance faith_pd_vector.log"
qiime diversity alpha-group-significance \
    --i-alpha-diversity "phylogenetic_core_metrics/evenness_vector.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --o-visualization "visualizations/alpha_evenness_group_significance.qzv" \
    --verbose \
    |& tee "${LOG_DIR}diversity alpha-group-significance evenness_vector.log"
qiime diversity alpha-rarefaction \
    --m-metadata-file "${METADATA_TSV}" \
    --i-table "dada2/dada2_frequency_table.qza" \
    --i-phylogeny "rooted_trees/rooted_tree.qza" \
    --o-visualization "visualizations/alpha_rarefaction.qzv" \
    --p-max-depth ${DENOISED_SAMPLES} \
    --verbose \
    |& tee "${LOG_DIR}diversity alpha-rarefaction.log"

log Visualize beta diversity
qiime diversity beta-group-significance \
    --i-distance-matrix "phylogenetic_core_metrics/unweighted_unifrac_distance_matrix.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --m-metadata-column "SampleSource" \
    --o-visualization "visualizations/beta_unweighted_unifrac_SampleSource_significance.qzv" \
    --p-pairwise \
    --verbose \
    |& tee "${LOG_DIR}diversity beta-group-significance SampleSource.log"
qiime diversity beta-group-significance \
    --i-distance-matrix "phylogenetic_core_metrics/unweighted_unifrac_distance_matrix.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --m-metadata-column "SubjectID" \
    --o-visualization "visualizations/beta_unweighted_unifrac_SubjectID_significance.qzv" \
    --p-pairwise \
    --verbose \
    |& tee "${LOG_DIR}diversity beta-group-significance SubjectID.log"
qiime emperor plot \
    --i-pcoa "phylogenetic_core_metrics/unweighted_unifrac_pcoa_results.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --o-visualization "visualizations/unweighted-unifrac-emperor.qzv" \
    --verbose \
    |& tee "${LOG_DIR}emperor plot.log"


log Test differential abundance with ANCOM
mkdir -p "differential_abundances/"
qiime feature-table filter-samples \
    --i-table "dada2/dada2_frequency_table.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --o-filtered-table "differential_abundances/source_features_table.qza" \
    --verbose \
    |& tee "${LOG_DIR}feature-table filter-samples.log"
qiime composition add-pseudocount \
    --i-table "differential_abundances/source_features_table.qza" \
    --o-composition-table "differential_abundances/pseudocounted_features_table.qza" \
    --verbose \
    |& tee "${LOG_DIR}composition add-pseudocount.log"
qiime composition ancom \
    --i-table "differential_abundances/pseudocounted_features_table.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --m-metadata-column "SubjectID" \
    --o-visualization "visualizations/ancom_SubjectID.qzv" \
    --verbose \
    |& tee "${LOG_DIR}composition ancom.log"


export COLLAPSE_LEVEL=6
log Test differential abundance with ANCOM for collapsed features
mkdir -p "collapsed_differential_abundances/"
qiime taxa collapse \
    --i-table "differential_abundances/source_features_table.qza" \
    --i-taxonomy "${TAXA_REFERENCE_FEATURES}" \
    --p-level ${COLLAPSE_LEVEL} \
    --o-collapsed-table "collapsed_differential_abundances/collapsed_source_features_table.qza" \
    --verbose \
    |& tee "${LOG_DIR}taxa collapse.log"
qiime composition add-pseudocount \
    --i-table "collapsed_differential_abundances/collapsed_source_features_table.qza" \
    --o-composition-table "collapsed_differential_abundances/collapsed_pseudocounted_features_table.qza" \
    --verbose \
    |& tee "${LOG_DIR}composition add-pseudocount collapsed.log"
qiime composition ancom \
    --i-table "collapsed_differential_abundances/collapsed_pseudocounted_features_table.qza" \
    --m-metadata-file "${METADATA_TSV}" \
    --m-metadata-column "SubjectID" \
    --o-visualization "visualizations/collapsed_ancom_SubjectID.qzv" \
    |& tee "${LOG_DIR}composition ancom collapsed.log"

log "Completed running QIIME2 in ${QIIME2_DIR}"
chmod -R 777 "$(pwd)"
cd ..
exit 0
