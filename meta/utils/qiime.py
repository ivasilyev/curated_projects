
import pandas as pd
from typing import Callable, List
from meta.sample_data.sample_data import create_sampledata_dict_from_dir


CATEGORICAL_TYPE = "categorical"
NUMERIC_TYPE = "numeric"
OTU_COLUMN_NAME = "#OTU ID"
SAMPLE_ID_NAME="#SampleID"
SAMPLE_ID_REGEX="[^A-Za-z0-9]+"


def fix_sample_ids(df: pd.DataFrame):
    # Deblur cannot operate on sample IDs that contain underscores
    df["#SampleID"].replace(
        SAMPLE_ID_REGEX,
        "-",
        inplace=True,
        regex=True,
    )


def split_metadata(df: pd.DataFrame) -> List[pd.DataFrame]:
    from meta.utils.pandas import split_df
    return split_df(df, 1, 0)


def fix_metadata(df: pd.DataFrame, sorting_column_name: str):
    from meta.utils.pandas import sort_df_by_values_count
    header_df, values_df = split_metadata(df)

    # Aggregating values sort
    sorted_values_df = sort_df_by_values_count(values_df, sorting_column_name)

    # Change sample ID
    fix_sample_ids(sorted_values_df)
    fixed_df = pd.concat(
        [header_df, sorted_values_df],
        axis=0,
        ignore_index=True,
        sort=False,
    )
    return fixed_df


def annotate_df_with_q2_types_row(df: pd.DataFrame):
    from numpy import issubdtype, number
    annotation_dict = dict()
    for column_name in df.columns:
        if column_name == "#SampleID":
            annotation_dict[column_name] = "#q2:types"
        elif issubdtype(df[column_name], number):
            annotation_dict[column_name] = NUMERIC_TYPE
        else:
            annotation_dict[column_name] = CATEGORICAL_TYPE
    return pd.concat(
        [pd.DataFrame([annotation_dict]), df],
        axis=0,
        sort=False,
    )

def create_main_metadata_df(
        directory: str,
        barcode_sequence: str = "",
        linker_primer_sequence: str = "",
        sample_source_extraction_function: Callable = None,
        subject_id_extraction_function: Callable = None,
        group_extraction_function: Callable = None,
        subgroup_extraction_function: Callable = None
):
    sample_data_dict = create_sampledata_dict_from_dir(directory)
    sample_data_dicts = list()
    for sampledata_line in sample_data_dict.values():
        sample_name = sampledata_line["name"]
        for sampledata_reads_file, direction in zip(
            sorted(sampledata_line["reads"]),
            ["forward", "reverse"]
        ):
            sample_data_dict = {
                SAMPLE_ID_NAME: sample_name,
                "BarcodeSequence": barcode_sequence,
                "LinkerPrimerSequence": linker_primer_sequence,
                "SampleSource": "",
                "ReadsStrand": direction,
                "SamplePath": sampledata_reads_file,
                "SubjectID": "",
                "Group": "",
                "Subgroup": ""
            }
            for func, key in zip([
                sample_source_extraction_function,
                subject_id_extraction_function,
                group_extraction_function,
                subgroup_extraction_function
            ], [
                "SampleSource",
                "SubjectID",
                "Group",
                "Subgroup"
            ]):
                if func is not None:
                    sample_data_dict[key] = func(sample_name)
            sample_data_dicts.append(sample_data_dict)
    raw_main_metadata_df = pd.DataFrame(sample_data_dicts).sort_values(SAMPLE_ID_NAME)
    main_metadata_df = annotate_df_with_q2_types_row(raw_main_metadata_df)
    return main_metadata_df


def convert_sampledata(
    sample_data_dict: dict,
    barcode_sequence: str = "",
    linker_primer_sequence: str = "",
):
    sample_data_dicts = []
    for sampledata_line in sample_data_dict.values():
        for sampledata_reads_file, direction in zip(
            sorted(sampledata_line["reads"]), ["forward", "reverse"]
        ):
            sample_data_dicts.append({
               "sample-id": sampledata_line["name"],
                "absolute-filepath": sampledata_reads_file,
                "direction": direction
            })
    meta_data_dicts = [{
        "#SampleID": "#q2:types",
        "BarcodeSequence": CATEGORICAL_TYPE,
        "LinkerPrimerSequence": CATEGORICAL_TYPE,
        "Description": CATEGORICAL_TYPE,
        "SampleSource": CATEGORICAL_TYPE,
        "SubjectID": CATEGORICAL_TYPE
    }, ]
    for sample_name in sorted(sample_data_dict.keys()):
        meta_data_dicts.extend([{
            "#SampleID": sample_name,
            "BarcodeSequence": barcode_sequence,
            "LinkerPrimerSequence": linker_primer_sequence,
            "Description": sample_name,
            "SampleSource": "".join([i for i in sample_name if i.isalpha()]),
            "SubjectID": "".join([i for i in sample_name if not i.isalpha()])
        }])
    return {
        "sample": pd.DataFrame(sample_data_dicts).sort_values("absolute-filepath"),
        "meta": pd.DataFrame(meta_data_dicts)
    }


def split_main_metadata_df(df: pd.DataFrame):
    main_metadata_sampledata_columns_mapper_dict = {
        SAMPLE_ID_NAME: "sample-id",
        "SamplePath": "absolute-filepath",
        "ReadsStrand": "direction",
    }
    sampledata_df = df.rename(
        main_metadata_sampledata_columns_mapper_dict,
        axis=1
    ).loc[
        :,
        main_metadata_sampledata_columns_mapper_dict.values()
    ].sort_values(
        main_metadata_sampledata_columns_mapper_dict[SAMPLE_ID_NAME]
    )
    _ = main_metadata_sampledata_columns_mapper_dict.pop(SAMPLE_ID_NAME)
    raw_metadata_df = df.drop(
        main_metadata_sampledata_columns_mapper_dict.keys(),
        axis=1
    ).drop_duplicates().sort_values(SAMPLE_ID_NAME)
    metadata_df = annotate_df_with_q2_types_row(raw_metadata_df)
    return dict(
        sample=sampledata_df,
        meta=metadata_df
    )


def dump_sampledata(sample_meta_data_dict: dict, directory: str):
    import os
    os.makedirs(directory, exist_ok=True)
    files_dict = dict()
    for key, df in sample_meta_data_dict.items():
        if key == "sample":
            sep = ","
            ext = "csv"
        else:
            sep = "\t"
            ext = "tsv"
        file = os.path.join(directory, f"qiime2_{key}_data.{ext}")
        df.to_csv(
            file,
            sep=sep,
            header=True,
            index=False
        )
        files_dict[key] = file
        print(f"Created {key} data: '{file}'")
    print(f"Sampledata created in: '{directory}'")
    return files_dict


def convert_and_dump_sampledata(directory: str, *args, **kwargs):
    dfs = convert_sampledata(*args, **kwargs)
    files_dict = dump_sampledata(dfs, directory)
    return files_dict
