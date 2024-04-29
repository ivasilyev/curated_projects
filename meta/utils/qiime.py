
import pandas as pd
from typing import List


_Q2_CAT_TYPE = "categorical"


def fix_sample_ids(df: pd.DataFrame):
    # Deblur cannot operate on sample IDs that contain underscores
    df["#SampleID"].replace(
        "[^A-Za-z0-9]+",
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
        "BarcodeSequence": _Q2_CAT_TYPE,
        "LinkerPrimerSequence": _Q2_CAT_TYPE,
        "Description": _Q2_CAT_TYPE,
        "SampleSource": _Q2_CAT_TYPE,
        "SubjectID": _Q2_CAT_TYPE
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


def convert_and_dump_sampledata(directory: str, *args, **kwargs):
    import os
    dfs = convert_sampledata(*args, **kwargs)
    os.makedirs(directory, exist_ok=True)
    out = dict()
    for key, df in dfs.items():
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
        out[key] = file
    print(f"Sampledata created in: '{directory}'")
    return out
