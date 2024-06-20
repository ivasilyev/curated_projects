#%%

import os
from meta.utils.pandas import load_tsv, dump_tsv
from meta.utils.file_system import filename_only, find_by_regex, get_file_extension


pipeline_output_dir = "/data03/bio/projects/ashestopalov/nutrition/mice_fatty_pentylresorcinol/qiime2-picrust2-pipeline"
kegg_reference_file = "/data/reference/KEGG/kegg_v2024-05-11/kegg_v2024-05-11_denormalized.tsv"
ko_file_regex = ".+qiime2.+dada2.+ASV.tsv$"

ko_sample_file = find_by_regex(".+picrust2.+KO_pred_metagenome_unstrat_described.tsv$", pipeline_output_dir)[0]

kegg_reference_df = load_tsv(kegg_reference_file)
kegg_reference_df

#%%

ko_sample_df = load_tsv(ko_sample_file)
ko_sample_df

#%%

ko_merged_df = ko_sample_df.merge(kegg_reference_df, how="left", on="function")
ko_merged_file_basename = f"{filename_only(ko_sample_file)}_annotated{get_file_extension(ko_sample_file)}"
ko_merged_file = os.path.join(os.path.dirname(ko_sample_file), ko_merged_file_basename)
dump_tsv(ko_merged_df, ko_merged_file)
