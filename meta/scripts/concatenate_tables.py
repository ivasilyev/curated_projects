#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from time import perf_counter
from meta.utils.pandas import load_tsv, dump_tsv
from meta.utils.file_system import is_file_valid
from meta.utils.primitive import remove_empty_values
from meta.utils.date_time import count_elapsed_seconds
from argparse import ArgumentParser, RawTextHelpFormatter


def parse_args():
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description="Concatenate tabular separated data which has same (or almost same) header",
        epilog=""
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
    tables = Utilities.remove_empty_values([i for i in input_tables if Utilities.is_file_valid(i)])
    if len(tables) == 0:
        raise ValueError("No valid tables!")
    dataframes = [Utilities.load_tsv(i) for i in tables]

    if len(index) > 0:
        for dataframe in dataframes:
            dataframe.set_index(index, inplace=True)

    print(f"Concatenate {len(dataframes)} dataframes with the shapes: {[i.shape for i in dataframes]}")
    start = perf_counter()
    out_df = pd.concat(dataframes, axis=axis, join=join, ignore_index=False, keys=None, levels=None,
                       names=None, verify_integrity=False, sort=False, copy=True).rename_axis(index=index)
    print(f"Concatenation completed in {Utilities.count_elapsed_seconds(start)}")

    Utilities.dump_tsv(out_df, output_table, reset_index=True)
