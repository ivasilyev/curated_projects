#%%

import os
import pandas as pd
from numpy import nan
from re import findall, sub
from meta.utils.file_system import find_by_regex, symlink
from meta.utils.language import translate_words
from meta.utils.pandas import (
    dfs_dict_to_excel,
    dump_tsv,
    excel_to_dfs_dict,
    split_df,
    set_df_dtypes)
from meta.utils.qiime import (
    create_main_metadata_df,
    dump_sampledata,
    SAMPLE_ID_NAME,
    split_main_metadata_df,
)

#%%

experimental_design_dfs_dict = excel_to_dfs_dict("/data03/bio/projects/ashestopalov/nutrition/mouse_obesity/sample_data/experimental_design.xlsx")
experimental_design_dfs_dict

#%%

group_df = experimental_design_dfs_dict["groups"]
group_df["GroupCode"] = group_df["group_code"].replace(
    "\-14$", "", regex=True
)
group_df.set_index("GroupCode", inplace=True)

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
        "group_code",
    ],
    axis=1,
    inplace=True
)
group_df

#%%

group_ids = set(group_df.index)
group_ids

#%%

source_raw_reads_dir = "/data03/bio/rogachev_human/blood/"
target_raw_reads_dir = "/data03/bio/rogachev_mice/raw/"

control_regex = "K[I]+"
control_reads = find_by_regex(f".+{control_regex}.+\.fastq\.gz$", source_raw_reads_dir)
symlinking_dicts = [
    {
        "source": i,
        "destination": os.path.join(target_raw_reads_dir, os.path.basename(i))
    }
    for i in control_reads
]

symlinking_dicts

#%%

_ = [symlink(**i) for i in symlinking_dicts]
control_symlinks = [i["destination"] for i in symlinking_dicts]
control_symlinks

#%%

main_metadata_df = create_main_metadata_df(
    #directory="/data03/bio/rogachev_mice/raw_colon_and_small/",
    directory=target_raw_reads_dir,
    barcode_sequence="CTCTCTACTATCCTCT",
    linker_primer_sequence="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG",
    sample_source_extraction_function=lambda x: "stool",
    subgroup_extraction_function=(
        lambda x: "066-" + sub("\-[^\-]+$", "", x)
    ),
    subject_id_extraction_function=lambda x: "066-" + x,
)

qiime_types, metadata_df = split_df(main_metadata_df, 1, 0)

excel_dict = {"all_files": metadata_df}
metadata_df

#%%

def validate_subject_id(subject_id: str):
    # Controls ex.: KIb, KIIb
    if len(findall(control_regex, subject_id)) > 0:
        return subject_id
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
feed_duration_column_name = "FeedDurationsDays"
metadata_df[feed_duration_column_name] = metadata_df.loc[:, "SubjectID"].str.extract("\-([^\-]+)\-[^\-]+$").loc[:, 0]
metadata_df["GroupCode"] = metadata_df["#SampleID"].map(
    lambda x: "066-" + sub("\-[^\-]+\-[^\-]+$", "", x)
)
metadata_df

#%%

merged_metadata_df = metadata_df.merge(
    group_df.reset_index(),
    how="left",
    left_on="GroupCode",
    right_on=group_df.index.name
).sort_values(SAMPLE_ID_NAME)

merged_metadata_df

#%%

sample_control_metadata_mapper = merged_metadata_df["SamplePath"].isin(control_symlinks)
control_metadata_df = merged_metadata_df.loc[
    sample_control_metadata_mapper,
    :
]
sample_metadata_df = merged_metadata_df.loc[
    ~sample_control_metadata_mapper,
    :
]
control_value = "ControlNegative"
for control_metadata_column_name in [
    "SampleSource",
    "FeedRation",
    "ExperimentCode",
    "Compound",
    "Subgroup",
    "GroupCode",
]:
   control_metadata_df.loc[:, control_metadata_column_name] = control_value

for control_metadata_column_name in [
    "Group",
    feed_duration_column_name
]:
   control_metadata_df.loc[:, control_metadata_column_name] = 0

control_metadata_df.loc[:, "SubjectID"] = control_metadata_df.loc[:, SAMPLE_ID_NAME]
control_metadata_df

#%%

sample_control_metadata_df = pd.concat(
    [sample_metadata_df, control_metadata_df,],
    axis=0,
    sort=True
)
set_df_dtypes(
    sample_control_metadata_df,
    {i: int for i in ["Group", feed_duration_column_name]}
)
sample_control_metadata_df

#%%

excel_dict.update({"Samplesheet": sample_control_metadata_df})
sampledata_dir = "/data03/bio/projects/ashestopalov/nutrition/mice_fatty_pentylresorcinol/sample_data/"
dfs_dict_to_excel(excel_dict, os.path.join(sampledata_dir, "samplesheet.xlsx"))
dump_tsv(sample_control_metadata_df, os.path.join(sampledata_dir, "main_sampledata.tsv"))

#%%

split_sampledata_dict = split_main_metadata_df(sample_control_metadata_df)
dump_sampledata(split_sampledata_dict, sampledata_dir)
split_sampledata_dict

#%%

# !head /data03/bio/projects/ashestopalov/nutrition/mice_fatty_pentylresorcinol/sample_data/qiime2_meta_data.tsv
# !tail /data03/bio/projects/ashestopalov/nutrition/mice_fatty_pentylresorcinol/sample_data/qiime2_meta_data.tsv

# Then upload the `main_sampledata.tsv`
