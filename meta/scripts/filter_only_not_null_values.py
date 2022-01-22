#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from numbers import Number


def is_value_not_null(x):
    if pd.isna(x):
        return False
    if isinstance(x, Number):
        return abs(x) > 0.0
    else:
        return len(str(x)) > 0


def parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description="Tool to keep only rows with the values from a column which are not null",
    )
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="Tables to concatenate")
    parser.add_argument("-f", "--filter", required=True,
                        help="Column name to filter")
    parser.add_argument("-o", "--output", required=True,
                        help="Output table")
    _namespace = parser.parse_args()
    return _namespace.input, _namespace.filter, _namespace.output


if __name__ == '__main__':
    input_file, filtering_column_name, output_file = parse_args()

    input_df = pd.read_csv(input_file, header=0, engine="python", on_bad_lines="warn", sep="\t")
    output_df = input_df.loc[input_df[filtering_column_name].map(is_value_not_null), :].sort_values(
        by=filtering_column_name, ascending=False)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    output_df.to_csv(output_file, encoding="utf-8", sep="\t", index=False, header=True)
