#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd


def load_tsv(table, col_names: list = None):
    if col_names:
        return pd.read_csv(table, encoding="utf-8", sep="\t", header="infer", names=col_names)
    return pd.read_csv(table, encoding="utf-8", sep="\t", header=0)


def dump_tsv(df: pd.DataFrame, table_file: str, col_names: list = None, reset_index: bool = False):
    assert isinstance(df, pd.DataFrame)
    _df = df.copy()
    os.makedirs(os.path.dirname(table_file), exist_ok=True)
    if col_names is not None and len(col_names) > 0:
        _df = _df.loc[:, col_names]
    if reset_index:
        _df.reset_index(inplace=True)
    _df.to_csv(table_file, encoding="utf-8", sep="\t", index=False, header=True)


def dict2pd_series(dictionary, sort_keys: bool = False):
    out = pd.Series()
    keys = list(dictionary.keys())
    if sort_keys:
        keys = sorted(keys)
    for key in keys:
        out.at[key] = dictionary[key]
    return out


def concat(dfs: list, index_name: str = "", columns_name: str = ""):
    return pd.concat(dfs, join="outer", axis=1, sort=False).rename_axis(
        index=index_name, columns=columns_name)
