#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from time import perf_counter
from meta.utils.file_system import is_file_valid
from meta.utils.primitive import remove_empty_values
from meta.utils.date_time import count_elapsed_seconds
from meta.utils.pandas import load_tsv, dump_tsv, dwell_df_on_column

PATH_COL_NAME = "source_table_file_name"


def parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description="Concatenate tabular separated data which has same (or almost same) header",
        epilog="The core method reference: https://pandas.pydata.org/docs/reference/api/pandas.concat.html"
    )
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="Tables to concatenate")
    parser.add_argument("-a", "--axis", default=0, type=int, choices=[0, 1],
                        help="(Optional) The axis to concatenate along")
    parser.add_argument("--index", default="",
                        help="(Optional) Index to concatenate with")
    parser.add_argument("--join", default="outer", choices=["inner", "outer"],
                        help="(Optional) How to handle indexes on other axis (or axes)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output table")
    _namespace = parser.parse_args()
    return _namespace.input, _namespace.axis, _namespace.index, _namespace.join, _namespace.output


if __name__ == '__main__':
    input_tables, axis, index, join, output_table = parse_args()
    table_files = remove_empty_values([i for i in input_tables if is_file_valid(i)])
    if len(table_files) == 0:
        raise ValueError("No valid tables!")
    dataframes = []
    for table_file in table_files:
        dataframe = load_tsv(table_file)
        if dataframe.shape[0] == 0:
            continue
        dataframe[PATH_COL_NAME] = table_file
        dataframes.append(dataframe)

    is_index = len(index) > 0
    if is_index:
        for dataframe in dataframes:
            dataframe.set_index(index, inplace=True)

    print(f"Concatenate {len(dataframes)} dataframes with the shapes: {[i.shape for i in dataframes]}")
    start = perf_counter()
    out_df = pd.concat(
        dataframes, axis=axis, join=join, ignore_index=False, keys=None, levels=None, names=None,
        verify_integrity=False, sort=False, copy=True
    ).rename_axis(index=index).sort_index()
    print(f"Concatenation completed in {count_elapsed_seconds(start)}")

    out_df = dwell_df_on_column(out_df, PATH_COL_NAME)
    dump_tsv(out_df, output_table, reset_index=is_index)
