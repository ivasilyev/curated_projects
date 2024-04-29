#%%

import os
import pandas as pd
from meta.utils.pandas import load_tsv, dump_tsv
from meta.utils.io import dump_list
from meta.utils.qiime import split_metadata


SAMPLE_DATA_COLUMN_DICT = {
    k: v
    for k, v in zip(
        "#SampleID	SamplePath	ReadsStrand".split("	"),
        "sample-id,absolute-filepath,direction".split(",")
    )
}
GROUP_COLUMN_NAME = "Subgroup"
NEGATIVE_CONTROL_GROUP_NAME = "ControlNegative"


#%%

def split_metadata_by_sample_group(
        sample_group: str,
        metadata_df: pd.DataFrame,
        output_dir: str
):
    """
    :param sample_group:
    :param metadata_df: DataFrame instance containing the first extra row with type descriptions
    :param output_dir:
    :return:
    """
    print(f"Processing '{sample_group}'")

    # Remove `#q2:types` row
    metadata_sample_df = metadata_df.drop([0], axis=0)

    group_sample_df = metadata_sample_df.loc[
        metadata_sample_df["Subgroup"].isin([sample_group, NEGATIVE_CONTROL_GROUP_NAME]),
        :
    ]

    group_meta_data_sample_df = group_sample_df.loc[
        :,
        [
            i for i in group_sample_df.columns
            if i not in list(SAMPLE_DATA_COLUMN_DICT.keys())[1:]
        ]
    ]
    group_meta_data_df = pd.concat(
        [
            pd.DataFrame(metadata_df.loc[
                :,
                group_meta_data_sample_df.columns
            ].iloc[0]).transpose(),
            group_meta_data_sample_df,
        ],
        axis=0,
        sort=False
    ).drop_duplicates("#SampleID")
    group_meta_data_file = os.path.join(output_dir, f"qiime2_meta_data-{sample_group}.tsv")
    dump_tsv(group_meta_data_df, group_meta_data_file, reset_index=False)
    print(f"Saved: '{group_meta_data_file}'")

    group_sample_data_df = group_sample_df.loc[
        :,
        SAMPLE_DATA_COLUMN_DICT.keys()
    ]
    group_sample_data_df = group_sample_data_df.rename(columns=SAMPLE_DATA_COLUMN_DICT)
    group_sample_data_file = os.path.join(output_dir, f"qiime2_sample_data-{sample_group}.csv")
    print(f"Saved: '{group_sample_data_file}'")
    group_sample_data_df.to_csv(group_sample_data_file, sep=",", index=False)


def split_and_dump_metadata(main_metadata_file: str, output_dir: str):
    main_metadata_df = load_tsv(main_metadata_file)

    _, values_df = split_metadata(main_metadata_df)

    sample_groups = [
        i for i in values_df[GROUP_COLUMN_NAME].unique()
        if i != NEGATIVE_CONTROL_GROUP_NAME
    ]

    for sample_group in sample_groups:
        split_metadata_by_sample_group(sample_group, main_metadata_df, output_dir)

    sample_groups_file = os.path.join(output_dir, "sample_groups.txt")

    dump_list(list(sample_groups), sample_groups_file)

    print(f"Saved sample groups: '{sample_groups_file}'")

#%%

split_and_dump_metadata(
    "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/stool_to_blood_2/sampledata/main_meta_data.tsv",
    "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/split/sample_data/"
)

#%%

# ! ls -lha /data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/sample_data/split

#%%

# ! head /data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/sample_data/split/qiime2_meta_data-NucChildren.tsv
