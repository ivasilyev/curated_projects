#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from meta.scripts.Utilities import Utilities
from argparse import ArgumentParser, RawTextHelpFormatter


def parse_args():
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description="Concatenate tabular separated data which has same (or almost same) header",
        epilog=""
    )
    parser.add_argument("-i", "--input", nargs="+",
                        help="Tables to concatenate")
    parser.add_argument("-a", "--axis", default=0, type=int, choices=[0, 1],
                        help="The axis to concatenate along.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output table")
    _namespace = parser.parse_args()
    return _namespace.input, _namespace.axis, _namespace.output


if __name__ == '__main__':
    input_tables, axis, output_table = parse_args()
    tables = Utilities.remove_empty_values([i for i in input_tables if Utilities.is_file_valid(i)])
    if len(tables) == 0:
        raise ValueError("No valid tables!")
    dfs = [Utilities.load_tsv(i) for i in tables]
    out_df = pd.concat(dfs, axis=axis, join="outer", ignore_index=False, keys=None, levels=None,
                       names=None, verify_integrity=False, sort=False, copy=True)
    Utilities.dump_tsv(out_df, output_table)
