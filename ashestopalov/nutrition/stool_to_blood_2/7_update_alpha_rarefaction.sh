
docker run --net=host --rm -it -v /data1:/data1 -v /data03:/data03 -v /data04:/data04 "quay.io/qiime2/core:2023.9"

# conda info --envs

source activate qiime2-amplicon-2023.9

conda config --set ssl_verify false
conda update \
    --channel defaults \
        conda \
    --name base \
    --yes \
    --verbose

SSL_NO_VERIFY=1 \
conda install \
    --channel bioconda \
    --channel conda-forge \
        bioconductor-phyloseq \
        r-devtools \
        r-tibble \
        r-magrittr \
        r-dplyr \
        r-withr \
        r-testthat \
        r-praise \
        unzip \
    --yes \
    --verbose


pip install git+https://github.com/statdivlab/q2-breakaway.git
R -e 'library(devtools); devtools::install_github("adw96/breakaway")'
qiime dev refresh-cache

qiime breakaway alpha --help

qiime breakaway alpha \
    --i-table /data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/denoising/frequency_table.qza \
    --o-alpha-diversity /data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/phylogenetic_core_metrics/richness-better.qza \
    --verbose

mkdir -p "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/phylogenetic_core_metrics/richness"

qiime tools export \
    --input-path "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/phylogenetic_core_metrics/richness-better.qza" \
    --output-path "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/phylogenetic_core_metrics/richness"

qiime breakaway plot \
    --i-alpha-diversity "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/phylogenetic_core_metrics/richness-better.qza" \
    --o-visualization "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/phylogenetic_core_metrics/richness/richness-better-plot.qzv" \
    --verbose

#    --i-table "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/vsearch/decontam/decontamination_filtered_frequencies.qza" \
mkdir -p "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/alpha_rarefaction"

qiime diversity alpha-rarefaction \
    --i-phylogeny "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/rooted_trees/rooted_tree.qza" \
    --i-table "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/denoising/frequency_table.qza" \
    --m-metadata-file "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/sample_data/qiime2_meta_data.tsv" \
    --o-visualization "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/alpha_rarefaction/alpha_rarefaction.gzv" \
    --p-max-depth 10501 \
    --p-metrics \
        berger_parker_d \
        brillouin_d \
        chao1 \
        dominance \
        doubles \
        enspie \
        faith_pd \
        fisher_alpha \
        gini_index \
        goods_coverage \
        heip_e \
        lladser_pe \
        margalef \
        mcintosh_d \
        mcintosh_e \
        menhinick \
        michaelis_menten_fit \
        observed_features \
        pielou_e \
        robbins \
        shannon \
        simpson \
        simpson_e \
        singles \
    --verbose
