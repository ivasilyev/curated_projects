
from meta.utils.pandas import dfs_dict_to_excel, load_tsv
from meta.utils.diversity import count_alpha_diversity_df
asv_table = "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/dada2/bioms/ASV_with_taxa.tsv"
asv_df = load_tsv(asv_table)
asv_df


#%%

indexed_asv_df = asv_df.drop(["taxonomy", "confidence"], axis=1).set_index("#OTU ID")
indexed_asv_df.columns.name = "#OTU ID"
asv_alpha_diversity_df = count_alpha_diversity_df(indexed_asv_df)
asv_alpha_diversity_df

#%%

otu_table = "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/qiime2/vsearch/bioms/OTU_normalized_with_taxa.tsv"
otu_df = load_tsv(otu_table)
otu_df

#%%

indexed_otu_df = otu_df.drop(["taxonomy",], axis=1).set_index("#OTU ID")
indexed_otu_df.columns.name = "#OTU ID"
otu_alpha_diversity_df = count_alpha_diversity_df(indexed_otu_df)
otu_alpha_diversity_df

#%%

alpha_diversity_table = "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/results/alpha_diversity.xlsx"
dfs_dict_to_excel(
    {
        "From ASV": asv_alpha_diversity_df.reset_index(),
        "From OTU": otu_alpha_diversity_df.reset_index()
    },
    alpha_diversity_table
)
