#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import pandas as pd
from meta.utils.primitive import remove_empty_values
from meta.utils.pandas import dfs_dict_to_excel, load_tsv, remove_longest_columns
from meta.utils.file_system import find_file_by_tail, filename_only, scan_top_level_directories


ANNOTATION_FILE_TAIL = "_coverage_annotated_filtered.tsv"
CELL_SIZE_LIMIT = 300


def _parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="")
    parser.add_argument("-r", "--rgi_dir", metavar="<dir>", default="",
                        help="RGI stage directory")
    parser.add_argument("-c", "--card_version", metavar="<str>", default="UNKNOWN",
                        help="CARD reference version")
    parser.add_argument("-n", "--nbee_dir", metavar="<dir>", required=True,
                        help="nBee stage directory")
    parser.add_argument("-o", "--output_file", metavar="<dir>", required=True,
                        help="Output file")
    _namespace = parser.parse_args()
    return (
        _namespace.rgi_dir,
        _namespace.card_version,
        _namespace.nbee_dir,
        _namespace.output_file
    )


if __name__ == '__main__':
    (
        rgi_dir,
        card_version,
        nbee_dir,
        out_file
    ) = _parse_args()
    sheets = dict()
    if len(rgi_dir) > 0:
        print("Process RGI")
        rgi_tables = find_file_by_tail(dir_name=rgi_dir, multiple=True, tail=".txt")
        merged_rgi_df = pd.DataFrame()
        for rgi_table in rgi_tables:
            rgi_df = load_tsv(rgi_table)
            if rgi_df.shape[0] == 0:
                continue
            rgi_df = remove_longest_columns(rgi_df, CELL_SIZE_LIMIT)
            columns = rgi_df.columns.tolist()
            rgi_df["sample_name"] = filename_only(rgi_table)
            rgi_df = rgi_df.loc[:, ["sample_name"] + columns]
            print(f"Concatenate dataframes with shapes {rgi_df.shape}, {merged_rgi_df.shape}")
            merged_rgi_df = pd.concat([merged_rgi_df, rgi_df], axis=0, ignore_index=True)
        reference_name = "_".join(remove_empty_values(["card", card_version]))
        print(f"Finished concatenating tables for '{reference_name}'")
        sheets[reference_name] = merged_rgi_df
    # Other references
    print("Process references")
    reference_dirs = scan_top_level_directories(nbee_dir)
    for reference_dir in reference_dirs:
        reference_name = os.path.basename(reference_dir)
        tail = f"_{reference_name}{ANNOTATION_FILE_TAIL}"
        coverage_tables = find_file_by_tail(
            dir_name=os.path.join(reference_dir, "annotated_coverages"), multiple=True, tail=tail
        )
        merged_coverage_df = pd.DataFrame()
        for coverage_table in coverage_tables:
            coverage_df = load_tsv(coverage_table)
            if coverage_df.shape[0] == 0:
                continue
            coverage_df = remove_longest_columns(coverage_df, CELL_SIZE_LIMIT)
            columns = coverage_df.columns.tolist()
            coverage_df["sample_name"] = os.path.basename(coverage_table).replace(tail, "")
            coverage_df = coverage_df.loc[:, ["sample_name"] + columns]
            print(f"Concatenate dataframes with shapes {coverage_df.shape}, {merged_coverage_df.shape}")
            merged_coverage_df = pd.concat([merged_coverage_df, coverage_df], axis=0, ignore_index=True)
        print(f"Finished concatenating tables for '{reference_name}'")
        sheets[reference_name] = merged_coverage_df
    dfs_dict_to_excel(sheets, out_file)
