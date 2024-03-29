#%%

from matplotlib import pyplot as plt
import os
import pandas as pd
import seaborn as sns
from time import perf_counter

from meta.sample_data.qiime2_sample_data import run
from meta.utils.file_system import find_file_by_tail
from meta.utils.file_system import scan_whole_dir
from meta.utils.diversity import (
    collapse_split_and_draw_non_major_df,
    count_alpha_diversity_df,
    count_beta_diversity_df,
    draw_dendrogram,
)
from meta.utils.pandas import (
    PRINCIPAL_ANALYSIS_METHODS,
    annotate_and_aggregate_df,
    count_feature_based_group_relations,
    dump_tsv,
    load_tsv,
)


def print_brief_df_report(df: pd.DataFrame):
    print("Columns:", df.columns.values, "\nShape:", df.shape)


start_timestamp = perf_counter()

# Raw folder is read-only!
raw_dir = "/data03/bio/rogachev_mice/raw_colon_and_small/"

# Output folders
project_dir =  "/data03/bio/projects/ashestopalov/nutrition/mice_small_large_intestine/"
qiime2_sample_data_dir = os.path.join(project_dir, "sample_data")
q2p2_pipeline_output_dir = os.path.join(project_dir, "results")
processed_data_dir = os.path.join(project_dir, "post-processed")
diversity_data_dir = os.path.join(processed_data_dir, "diversity")
pivot_diversity_data_dir = os.path.join(diversity_data_dir, "pivots")
dendrogram_dir = os.path.join(diversity_data_dir, "dendrograms")
collapsed_diversity_data_dir = os.path.join(diversity_data_dir, "collapsed-taxa")
aggregated_diversity_data_dir = os.path.join(diversity_data_dir, "aggregated-taxa")
pathway_data_dir = os.path.join(processed_data_dir, "pathways")
aggregated_pathway_data_dir = os.path.join(pathway_data_dir, "aggregated-pathways")

#%%

raw_reads_files = pd.Series(scan_whole_dir(raw_dir))
raw_reads_files

#%%

qiime2_run_dicts = run(
    reads_dirs=[raw_dir,],
    extension=".fastq.gz",
    sampledata_dir=qiime2_sample_data_dir,
    # Illumina_16S_341F primer for V3-V4 16S region
    barcode_sequence="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG",
    linker_primer_sequence="",
)

qiime2_metadata_file = qiime2_run_dicts["sampledata_files"]["meta"]
qiime2_metadata_df = load_tsv(qiime2_metadata_file).iloc[1:]

#%%

index_column_name = "#SampleID"
qiime2_metadata_df["SampleSource"] = qiime2_metadata_df[index_column_name].map(
    lambda x:
    "large_intestine" if any(x.endswith(i) for i in ["c", "ls"])
    else "small_intestine"

)

sample_names = qiime2_metadata_df[index_column_name]
qiime2_metadata_df["SubjectID"] = sample_names.str.extract(
    "([0-9]+\-[0-9]+)", expand=False
).values.tolist()

#%%

qiime2_metadata_df = pd.concat([
    load_tsv(qiime2_metadata_file).iloc[:1],
    qiime2_metadata_df
], axis=0)

dump_tsv(qiime2_metadata_df, qiime2_metadata_file)

#%%

qiime2_sampledata_df = pd.read_csv(qiime2_run_dicts["sampledata_files"]["sample"])
qiime2_sampledata_df["size_bytes"] = qiime2_sampledata_df["absolute-filepath"].apply(
    lambda x: os.stat(x).st_size
)
qiime2_sampledata_df

#%%

plt.clf()
plt.close()

sns.set()

plt.rcParams.update({
    "figure.figsize": (5, 5),
    "figure.dpi": 75
})

ax = sns.violinplot(
    data=qiime2_sampledata_df,
    x="direction",
    y="size_bytes",
    palette="YlOrBr"
)
plt.suptitle("Raw file size distribution")
# ax.set(ylabel="size_bytes")
plt.show()

#%%

_ = """
export ROOT_DIR="/data03/bio/projects/ashestopalov/nutrition/mice_small_large_intestine/"
export RAW_DIR="/data03/bio/rogachev_mice/raw_colon_and_small/"

export SAMPLEDATA_DIR="${ROOT_DIR}sample_data/"
export SCRIPT_DIR="${ROOT_DIR}scripts/"
export PIPELINE_SCRIPT="${SCRIPT_DIR}1_run_pipeline"

# rm -rf "${ROOT_DIR}"
mkdir -p "${ROOT_DIR}" "${SCRIPT_DIR}"
chmod -R a+rw "${ROOT_DIR}"
cd "${ROOT_DIR}" || exit 1

curl -fsSL "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/mouse_obesity/1_run_pipeline.sh" \
    -o "${PIPELINE_SCRIPT}"

ROOT_DIR="${ROOT_DIR}" \
RAW_DIR="${RAW_DIR}" \
SAMPLEDATA_DIR="${SAMPLEDATA_DIR}" \
bash "${PIPELINE_SCRIPT}"
"""

#%%

otu_df = load_tsv(
    find_file_by_tail(q2p2_pipeline_output_dir, "OTUs_with_taxa_annotated.tsv")
).set_index("#OTU ID").rename_axis(columns=index_column_name)
sample_otu_df = otu_df.loc[:, sample_names].fillna(0)


#%%

alpha_diversity_df = count_alpha_diversity_df(sample_otu_df)
dump_tsv(
    alpha_diversity_df,
    os.path.join(pivot_diversity_data_dir, "alpha_diversity_table.tsv"),
    reset_index=True
)
alpha_diversity_df

#%%

relation_feature_name = "OTU"

source_column_name = "SampleSource"
source_otu_grouping_df = pd.concat(
    [
        qiime2_metadata_df.iloc[1:].loc[
            :,
            [index_column_name, source_column_name]
        ].set_index(index_column_name),
        sample_otu_df.transpose()
    ],
    axis=1,
    sort=False,
    join="inner"
)

#%%

os.environ["MATPLOTLIB_COLORMAP"] = "hsv_r"
source_otu_count_df_dict = count_feature_based_group_relations(
    df=source_otu_grouping_df,
    grouping_column_name=source_column_name,
    feature_name=relation_feature_name,
    output_dir=diversity_data_dir,
    annotation_df=pd.DataFrame(otu_df["taxonomy"])
)

#%%

os.environ["MATPLOTLIB_COLORMAP"] = "turbo_r"

for principal_method_name, principal_method_function in PRINCIPAL_ANALYSIS_METHODS.items():
    principal_method_function(
        df=source_otu_grouping_df,
        grouping_column_name=source_column_name,
        output_dir=os.path.join(
            diversity_data_dir,
            principal_method_name,
            "{}_for_{}_by_{}".format(
                principal_method_name,
                relation_feature_name,
                source_column_name
            )
        ),
        name=f"{relation_feature_name} grouped by {source_column_name}",
    )

#%%

pathway_df = load_tsv(
    find_file_by_tail(
        q2p2_pipeline_output_dir,
        "KO_metagenome_out_pred_metagenome_unstrat_described.tsv"
    ),
    index_col="function"
)

#%%

relation_feature_name = "enzyme"
os.environ["MATPLOTLIB_COLORMAP"] = "hsv"
source_enzyme_grouping_df = pd.concat(
    [
        qiime2_metadata_df.iloc[1:].loc[
            :,
            [index_column_name, source_column_name]
        ].set_index(index_column_name),
        pathway_df.loc[:, sample_names].transpose()
    ],
    axis=1,
    sort=False
).rename_axis(index=index_column_name, columns=relation_feature_name)

#%%

source_enzyme_count_df_dict = count_feature_based_group_relations(
    df=source_enzyme_grouping_df,
    grouping_column_name=source_column_name,
    feature_name=relation_feature_name,
    output_dir=pathway_data_dir,
    annotation_df=pd.DataFrame(pathway_df["description"])
)

#%%

os.environ["MATPLOTLIB_COLORMAP"] = "turbo"
for principal_method_name, principal_method_function in PRINCIPAL_ANALYSIS_METHODS.items():
    principal_method_function(
        df=source_enzyme_grouping_df,
        grouping_column_name=source_column_name,
        output_dir=os.path.join(
            pathway_data_dir,
            principal_method_name,
            "{}_for_{}_by_{}".format(
                principal_method_name,
                relation_feature_name,
                source_column_name
            )
        ),
        name=f"{relation_feature_name} grouped by {source_column_name}",
    )

#%%

pathway_reference_dir = "/data/reference/KEGG/"
pathway_reference_table = find_file_by_tail(pathway_reference_dir, ".tsv")
pathway_reference_df = load_tsv(pathway_reference_table).set_index("ko")
pathway_reference_df

#%%

pre_aggregated_pathway_df, aggregated_pathway_dfs_dict = annotate_and_aggregate_df(
    df=pathway_df.loc[:, sample_names].rename_axis(index="ko"),
    annotation_df=pathway_reference_df,
)

#%%

dump_tsv(
    pre_aggregated_pathway_df,
    os.path.join(aggregated_pathway_data_dir, "pre_aggregated_pathways.tsv"),
    reset_index=True,
)

_ = [
    dump_tsv(
        v,
        os.path.join(aggregated_pathway_data_dir, f"pathways_aggregated_by_{k}.tsv"),
        reset_index=True,
    )
    for k, v in aggregated_pathway_dfs_dict.items()
]

#%%

subject_column_name = "SubjectID"
subject_otu_grouping_df = pd.concat(
    [
        qiime2_metadata_df.iloc[1:].loc[
            :,
            [index_column_name, subject_column_name]
        ].set_index(index_column_name),
        sample_otu_df.transpose()
    ],
    axis=1,
    sort=False
)

#%%

beta_diversity_df = count_beta_diversity_df(
    df=subject_otu_grouping_df,
    grouping_column_name=subject_column_name
)
dump_tsv(
    beta_diversity_df,
    os.path.join(pivot_diversity_data_dir, "beta_diversity_table.tsv"),
    reset_index=True
)
beta_diversity_df

#%%

draw_dendrogram(
    df=sample_otu_df.loc[
        :,
        [
            i for i in sample_otu_df
            if all(not i.startswith(j) for j in ["7", "8",])
        ]
    ].transpose(),
    name="OTU-based",
    metric="braycurtis",
    output_dir=dendrogram_dir,
)

#%%

taxa_sample_otu_df = otu_df.loc[:, ["taxonomy", *sample_names]]
taxa_sample_otu_df

#%%

collapse_split_and_draw_non_major_df(
    taxa_sample_otu_df=taxa_sample_otu_df,
    major_features_number=20,
    samples_per_time=10,
    output_dir=collapsed_diversity_data_dir
)

#%%

taxa_reference_dir = "/data/reference/SILVA/SILVA_v138"
taxa_reference_table = find_file_by_tail(taxa_reference_dir, "_taxonomy_expanded.tsv")
taxa_reference_df = load_tsv(taxa_reference_table).set_index("#OTU ID").drop("taxonomy", axis=1)
taxa_reference_df

#%%

pre_aggregated_taxa_df, aggregated_taxa_dfs_dict = annotate_and_aggregate_df(
    df=otu_df.loc[:, sample_names],
    annotation_df=taxa_reference_df,
)

#%%

dump_tsv(
    pre_aggregated_taxa_df,
    os.path.join(aggregated_diversity_data_dir, "pre_aggregated_taxa.tsv"),
    reset_index=True,
)

_ = [
    dump_tsv(
        v,
        os.path.join(aggregated_diversity_data_dir, f"taxa_aggregated_by_{k}.tsv"),
        reset_index=True,
    )
    for k, v in aggregated_taxa_dfs_dict.items()
]

#%%

print(f"Time passed: {perf_counter() - start_timestamp:.3f} s.")
