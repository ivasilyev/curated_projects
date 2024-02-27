#%%

import os
import pandas as pd
from meta.utils.pandas import load_tsv, dump_tsv


SAMPLE_DATA_COLUMN_DICT = {
    k: v
    for k, v in zip(
        "#SampleID	SamplePath	ReadsStrand".split("	"),
        "sample-id,absolute-filepath,direction".split(",")
    )
}


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
    metadata_sample_df = metadata_df.drop([0], axis=0)

    group_sample_df = metadata_sample_df.loc[
        metadata_sample_df["SampleGroup1"] == sample_group,
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
            pd.DataFrame(metadata_sample_df.loc[
                :,
                group_meta_data_sample_df.columns
            ].iloc[0]).transpose(),
            group_meta_data_sample_df,
        ],
        axis=0,
        sort=False
    ).drop_duplicates("#SampleID")
    dump_tsv(
        group_meta_data_df,
        os.path.join(output_dir, f"qiime2_meta_data-{sample_group}.tsv"),
        reset_index=False
    )

    group_sample_data_df = group_sample_df.loc[
        :,
        SAMPLE_DATA_COLUMN_DICT.keys()
    ]
    group_sample_data_df = group_sample_data_df.rename(columns=SAMPLE_DATA_COLUMN_DICT)
    group_sample_data_df.to_csv(
        os.path.join(
            output_dir,
            f"qiime2_sample_data-{sample_group}.csv"
        ),
        sep=",",
        index=False
    )


def split_metadata(main_metadata_file: str, output_dir: str):
    main_metadata_df = load_tsv(main_metadata_file)

    for sample_group in set(main_metadata_df.drop([0], axis=0)["SampleGroup1"].values):
        split_metadata_by_sample_group(sample_group, main_metadata_df, output_dir)

#%%

split_metadata(
    "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/ashestopalov/nutrition/stool_to_blood_2/sampledata/main_meta_data.tsv",
    "/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/sample_data/split"
)

#%%

# ! ls /data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2/qiime2-picrust2-pipeline/sample_data/split
