#%%

import pandas as pd
from numpy import nan
from re import sub
from meta.utils.qiime import create_main_metadata_df
from meta.utils.pandas import dfs_dict_to_excel, excel_to_dfs_dict, split_df
from meta.utils.language import translate_words

#%%

experimental_design_dfs_dict = excel_to_dfs_dict("/data03/bio/projects/ashestopalov/nutrition/mouse_obesity/sample_data/experimental_design.xlsx")
experimental_design_dfs_dict.keys()

#%%

group_df = experimental_design_dfs_dict["groups"].set_index("group_code")

ration_column_name = "ration"
en_rations = translate_words(group_df[ration_column_name].values.tolist())

group_df[ration_column_name] = group_df.loc[:, ration_column_name].replace(en_rations)

group_df = group_df.rename(
    columns={i: i.strip() for i in group_df.columns}
).rename(
    columns={
        "experiental_subgroup": "Subgroup",
        "mouse_count": "SubjectsPerGroupCount",
        "experimental_code": "ExperimentCode",
        "ration": "FeedRation",
        "days_on_ration": "FeedDurationsDays",
    }
)

group_df["Compound"] = group_df["is_solvent"].map(
    lambda x: "solvent" if x == "+" else "pentylresorcinol"
)

group_df.drop(
    ["is_pentylresorcinol", "is_solvent"],
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
    sample_source_extraction_function= lambda x: "stool",
    group_extraction_function=lambda x: "066-" + sub("\-[^\-]+$", "", x),
    subject_id_extraction_function=lambda x: x,
)

qiime_types, metadata_df = split_df(main_metadata_df, 1, 0)

excel_dict = {"all_files": metadata_df}

metadata_df

#%%

def validate_sample_name(sample_name: str):
    if sample_name.count("-") != 2:
        return nan
    group_number, elapsed_days_number, sample_number = sample_name.split("-")
    full_group_id = f"066-{group_number}-{elapsed_days_number}"
    if full_group_id not in group_ids:
        return nan
    return f"066-{sample_name}"


metadata_df["#SampleID"] = metadata_df["#SampleID"].map(validate_sample_name)
metadata_df = metadata_df.dropna(axis=0, how="any").drop(["Subgroup"], axis=1)
metadata_df

#%%


merged_metadata_df = metadata_df.merge(
    group_df.reset_index(),
    how="left",
    left_on="Group",
    right_on=group_df.index.name
).drop(
    [group_df.index.name,],
    axis=1
)
merged_metadata_df

#%%

_Q2_CAT_TYPE = "categorical"


def annotate_df_with_qiime_row(df: pd.DataFrame):
    dd = [{
            i: "#q2:types"
            if i == "#SampleID"
            else _Q2_CAT_TYPE
            for i in df.columns
    }]
    return pd.concat(
        [pd.DataFrame(dd), df],
        axis=0,
        sort=False,
    )


annotated_merged_metadata_df = annotate_df_with_qiime_row(merged_metadata_df)
annotated_merged_metadata_df

#%%

excel_dict.update({"Samplesheet": annotated_merged_metadata_df})

dfs_dict_to_excel(excel_dict, "/data03/bio/projects/ashestopalov/nutrition/mouse_obesity/sample_data/samplesheet.xlsx")
