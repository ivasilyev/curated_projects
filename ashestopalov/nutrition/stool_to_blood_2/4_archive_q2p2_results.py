#%%

import os
import pandas as pd
from meta.utils.file_system import filename_only, scan_whole_dir
from meta.utils.diversity import annotate_df_with_taxa_columns
from meta.utils.pandas import annotate_and_aggregate_df, df_to_7z, load_tsv, split_df
from meta.utils.primitive import dicts_list_to_lists_dict
from meta.utils.queue import multi_core_queue2

#%%

asset_basenames = {
    "ASV.tsv",
    "ASV_confidences.tsv",
    "OTU_with_taxa_normalized.tsv",
    # "EC_pred_metagenome_contrib_legacy.tsv",  # OOMKilled
    "EC_pred_metagenome_unstrat_described.tsv",
    "KO_pred_metagenome_unstrat_described.tsv",
    "pathways_abun_unstrat_described.tsv",
}

pipeline_output_dir = "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline"

pipeline_output_files = {
    k: v
    for k in asset_basenames
    for v in scan_whole_dir(pipeline_output_dir)
    if v.endswith(k)
}
pipeline_output_files

#%%

pipeline_output_dfs_list = multi_core_queue2(
    lambda x: {os.path.basename(x): load_tsv(x)},
    pipeline_output_files.values()
)
pipeline_output_dfs_dict = dicts_list_to_lists_dict(pipeline_output_dfs_list)

#%%

ko_reference_table = "/data/reference/KEGG/kegg_v2024-05-11/kegg_v2024-05-11_denormalized.tsv"
ko_reference_df = load_tsv(ko_reference_table, index_col="function")
ko_reference_df

#%%


ko_raw_df = pipeline_output_dfs_dict["KO_pred_metagenome_unstrat_described.tsv"][0]
ko_annotated_df, ko_aggregated_dfs_dict = annotate_and_aggregate_df(
    df=ko_raw_df.set_index(["function", "description"]),
    annotation_df=ko_reference_df,
)

ko_annotated_df

#%%

pipeline_output_dfs_dict["KO_annotated.tsv"] = [ko_annotated_df.reset_index()]

#%%

asv_reference_df = pipeline_output_dfs_dict["ASV_confidences.tsv"][0].set_index("#OTU ID")
asv_reference_df

#%%

asv_reference_taxa_df = annotate_df_with_taxa_columns(asv_reference_df)
asv_reference_taxa_df

#%%

asv_confidence_df, asv_annotation_df = split_df(asv_reference_taxa_df, 2, 1)
asv_confidence_df

#%%

asv_raw_df = pipeline_output_dfs_dict["ASV.tsv"][0]

if "taxonomy" in asv_raw_df.columns:
    asv_raw_df.drop("taxonomy", axis=1, inplace=True)

asv_annotated_df, asv_aggregated_dfs_dict = annotate_and_aggregate_df(
    df=asv_raw_df.set_index("#OTU ID"),
    annotation_df=asv_annotation_df,
)
asv_annotated_df

#%%

asv_annotated_confidence_df = pd.concat([asv_confidence_df, asv_annotated_df], axis=1, sort=False)
asv_annotated_confidence_df

#%%

pipeline_output_dfs_dict["ASV_annotated.tsv"] = [asv_annotated_confidence_df.reset_index()]

#%%

results_dir = "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/results"

def export_assets(pipeline_output_dfs_dict_key):
    df = pipeline_output_dfs_dict.get(pipeline_output_dfs_dict_key)[0]
    table_basename = filename_only(pipeline_output_dfs_dict_key)
    archive = os.path.join(results_dir, f"{table_basename}.7z")
    print(f"Processing '{pipeline_output_dfs_dict_key}'")
    df_to_7z(df=df, archive=archive, table_basename=table_basename)


multi_core_queue2(
    export_assets,
    list(pipeline_output_dfs_dict.keys())
)
