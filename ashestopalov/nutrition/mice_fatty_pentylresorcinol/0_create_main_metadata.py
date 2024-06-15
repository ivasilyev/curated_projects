#%%

import os
from numpy import nan
from re import sub
from meta.utils.qiime import create_main_metadata_df, dump_sampledata, split_main_metadata_df
from meta.utils.pandas import dfs_dict_to_excel, dump_tsv, excel_to_dfs_dict, split_df
from meta.utils.language import translate_words

#%%

experimental_design_dfs_dict = excel_to_dfs_dict("/data03/bio/projects/ashestopalov/nutrition/mouse_obesity/sample_data/experimental_design.xlsx")
experimental_design_dfs_dict

#%%

group_df = experimental_design_dfs_dict["groups"]
group_df["group_code"] = group_df["group_code"].replace(
    "\-14$", "", regex=True
)
group_df.set_index("group_code", inplace=True)

ration_column_name = "ration"
en_rations = translate_words(group_df[ration_column_name].values.tolist())

group_df[ration_column_name] = group_df.loc[:, ration_column_name].replace(en_rations)

group_df = group_df.rename(
    columns={i: i.strip() for i in group_df.columns}
).rename(
    columns={
        "experiental_subgroup": "Group",
        "experimental_code": "ExperimentCode",
        "ration": "FeedRation",
    }
)

group_df["Compound"] = group_df["is_solvent"].map(
    lambda x: "solvent" if x == "+" else "pentylresorcinol"
)
group_df.drop(
    [
        "days_on_ration",
        "is_pentylresorcinol",
        "is_solvent",
        "mouse_count",
    ],
    axis=1,
    inplace=True
)
group_df

#%%

group_ids = set(group_df.index)
group_ids

#%%

main_metadata_df = create_main_metadata_df(
    #directory="/data03/bio/rogachev_mice/raw_colon_and_small/",
    directory="/data03/bio/rogachev_mice/raw/",
    barcode_sequence="CTCTCTACTATCCTCT",
    linker_primer_sequence="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG",
    sample_source_extraction_function=lambda x: "stool",
    subgroup_extraction_function=(
        lambda x: "066-" + sub("\-[^\-]+\-[^\-]+$", "", x)
    ),
    subject_id_extraction_function=lambda x: "066-" + x,
)

qiime_types, metadata_df = split_df(main_metadata_df, 1, 0)

excel_dict = {"all_files": metadata_df}

metadata_df

#%%

def validate_subject_id(subject_id: str):
    # Ex.: 066-050-14-10
    if subject_id.count("-") != 3:
        return nan
    prefix, group_number, elapsed_days_number, sample_number = subject_id.split("-")
    full_group_id = f"{prefix}-{group_number}"
    if full_group_id not in group_ids:
        return nan
    return subject_id


metadata_df["ValidSubjectID"] = metadata_df.loc[:, "SubjectID"].map(validate_subject_id)
metadata_df = metadata_df.dropna(axis=0, how="any").drop(["ValidSubjectID", "Group"], axis=1)
metadata_df["FeedDurationsDays"] = metadata_df.loc[:, "SubjectID"].str.extract("\-([^\-]+)\-[^\-]+$").loc[:, 0]
metadata_df

#%%

merged_metadata_df = metadata_df.merge(
    group_df.reset_index(),
    how="left",
    left_on="Subgroup",
    right_on=group_df.index.name
).drop(
    [group_df.index.name,],
    axis=1
)
merged_metadata_df["FeedDurationsDays"] = merged_metadata_df.loc[:, "FeedDurationsDays"].astype(int)
merged_metadata_df

#%%

excel_dict.update({"Samplesheet": merged_metadata_df})
sampledata_dir = "/data03/bio/projects/ashestopalov/nutrition/mice_fatty_pentylresorcinol/sample_data/"
dfs_dict_to_excel(excel_dict, os.path.join(sampledata_dir, "samplesheet.xlsx"))
dump_tsv(merged_metadata_df, os.path.join(sampledata_dir, "main_sampledata.tsv"))

#%%

split_sampledata_dict = split_main_metadata_df(merged_metadata_df)

#%%

dump_sampledata(split_sampledata_dict, sampledata_dir)
