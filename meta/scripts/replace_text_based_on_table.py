#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from meta.utils.pandas import load_tsv
from meta.utils.io import load_string, dump_string


def _parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description="Replace GenInfo ID with strains from given table".strip(),
    )
    parser.add_argument("--target_file", required=True,
                        help="Target file with text to perform replacements")
    parser.add_argument("--target_column", required=True,
                        help="Column name with values to be replaced")
    parser.add_argument("--source_table", required=True,
                        help="Source table with replacers")
    parser.add_argument("--source_columns", required=True, nargs="+",
                        help="Column name(s) containing replacer(s)")
    parser.add_argument("--separator", default=" ",
                        help="Separator to join, the space as default")
    parser.add_argument("o", "--output_file", required=True,
                        help="Output file")
    _namespace = parser.parse_args()
    return (
        _namespace.target_file,
        _namespace.target_column,
        _namespace.source_table,
        _namespace.source_columns,
        _namespace.separator,
        _namespace.output_file,
    )


if __name__ == '__main__':
    (
        target_file,
        target_column_name,
        source_table_file,
        source_column_names,
        separator,
        out_file,
    ) = _parse_args()

    source_df = load_tsv(source_table_file)
    source_sub_df = pd.DataFrame()
    source_sub_df["target"] = source_df[target_column_name].copy()
    source_sub_df["source"] = source_df.loc[:, source_column_names].fillna("").agg(separator.join, axis=1)
    source_sub_df = source_sub_df.loc[source_sub_df["source"].map(lambda x: len(x) > 0), :]

    target_string = load_string(target_file)

    collector_string = str(target_string)
    counter = 0
    for source, target in source_sub_df.values:
        collector_string = collector_string.replace(str(source), target)
        counter += 1

    dump_string(collector_string, out_file)
    print(f"Saved: '{out_file}' with {counter} replacements made")
