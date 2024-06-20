#%%

import os
from meta.utils.pandas import dump_tsv, load_tsv
from meta.utils.file_system import filename_only, find_by_regex, get_file_extension


pipeline_output_dir = "/data03/bio/projects/ashestopalov/nutrition/mice_fatty_pentylresorcinol/qiime2-picrust2-pipeline"
asv_file_regex = ".+qiime2.+dada2.+ASV.tsv$"
mapper_file_regex = ".+qiime2.+dada2.+ASV_confidences.tsv$"
index_column_name = "#OTU ID"

asv_raw_file = find_by_regex(asv_file_regex, pipeline_output_dir)[0]
asv_raw_df = load_tsv(asv_raw_file)
asv_df = asv_raw_df.drop(["taxonomy"], axis=1)
asv_df

#%%

otu_asv_mapper_file = find_by_regex(mapper_file_regex, pipeline_output_dir)[0]
otu_asv_mapper_df = load_tsv(otu_asv_mapper_file)
otu_asv_mapper_df

#%%

asv_annotated_df = otu_asv_mapper_df.merge(
    asv_df,
    how="right",
    on=index_column_name
)
# assert not asv_annotated_df[index_column_name].duplicated().any()
asv_annotated_df

#%%

asv_merged_file_basename = f"{filename_only(asv_raw_file)}_with_taxa{get_file_extension(asv_raw_file)}"
asv_merged_file = os.path.join(os.path.dirname(asv_raw_file), asv_merged_file_basename)
print(asv_merged_file)
dump_tsv(asv_annotated_df, asv_merged_file)
