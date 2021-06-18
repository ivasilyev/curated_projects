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


def apply_mp_function_to_df(func, df: pd.DataFrame, index_name: str = "", columns_name: str = ""):
    from meta.scripts.utils.queue_utils import multi_core_queue
    results = multi_core_queue(func, [df[i] for i in df.columns], async_=True)
    return pd.DataFrame(results).rename_axis(index=index_name, columns=columns_name)


def corr(df: pd.DataFrame, methods: list = None):
    from numpy import eye
    from scipy import stats
    if methods is None or len(methods) != 2:
        methods = [lambda x, y: stats.spearmanr(x, y)[0], lambda x, y: stats.spearmanr(x, y)[-1]]
    correlation_df = df.corr(method=methods[0])
    p_values_df = df.corr(method=methods[-1]) - eye(*correlation_df.shape)
    asterisk_df = p_values_df.applymap(
        lambda x: "".join(["*" for i in (0.01, 0.05, 0.1) if x < i]))
    denoted_correlation_df = correlation_df.round(2).astype(str) + asterisk_df
    d = {"correlations": correlation_df, "p_values": p_values_df,
         "denoted_correlation": denoted_correlation_df}
    try:
        d["name"] = df.name
    except AttributeError:
        pass
    return d
