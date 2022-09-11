#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import joblib
import traceback
import numpy as np
import pandas as pd
from scipy import stats


def load_tsv(table, **kwargs):
    _kwargs = dict(encoding="utf-8", engine="python", header=0, sep="\t")
    if len(kwargs.keys()) > 0:
        _kwargs.update(kwargs)
    return pd.read_csv(table, **_kwargs)


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


def apply_mp_function_to_df(func, df: pd.DataFrame, index_name: str = "", columns_name: str = ""):
    results = joblib.Parallel(n_jobs=-1)(joblib.delayed(func)(j) for j in [df[i] for i in df.columns])
    return pd.DataFrame(results).rename_axis(index=index_name, columns=columns_name)


def corr(df: pd.DataFrame, methods: list = None):
    if methods is None or len(methods) != 2:
        methods = [lambda x, y: stats.spearmanr(x, y)[0], lambda x, y: stats.spearmanr(x, y)[-1]]
    correlation_df = df.corr(method=methods[0])
    p_values_df = df.corr(method=methods[-1]) - np.eye(*correlation_df.shape)
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


def convert_objects(df: pd.DataFrame, type_=str):
    _df = df.copy()
    for column in _df.columns:
        _df[column] = _df[column].astype(type_)
    return _df


def find_duplicated_rows(df: pd.DataFrame, on: str):
    return df.loc[df[on].duplicated(keep=False)].sort_values(on)


def find_notna_rows(df: pd.DataFrame, on: str):
    return df.loc[~df[on].isna()].sort_values(on)


def deduplicate_df_by_row_merging(df: pd.DataFrame, on: str, sep: str = ";"):
    def _deduplicate(series: pd.Series):
        return sep.join(sorted(set([str(i) for i in series.values])))
    _df = convert_objects(df)
    duplicated_df = find_duplicated_rows(_df, on)
    out_series = []
    for duplicated_value in sorted(set(duplicated_df[on].values)):
        out_series.append(
            duplicated_df.loc[duplicated_df[on] == duplicated_value].apply(_deduplicate, axis=0)
        )
    return pd.concat(
        [pd.DataFrame(out_series), _df.drop(index=duplicated_df.index)],
        axis=0
    ).sort_values(on)


def concat(dfs: list, on: str, axis: int = 0, how: str = "outer", columns_name: str = "",
           set_index: bool = True, reset_index: bool = True, **kwargs):
    assert all(isinstance(i, pd.DataFrame) for i in dfs)
    _dfs = list(dfs)
    if set_index:
        _dfs = [i.set_index(on) for i in _dfs]
    out = pd.concat(
        _dfs, join=how, axis=axis, **kwargs
    ).rename_axis(index=on, columns=columns_name).sort_index()
    if reset_index:
        return out.reset_index()
    return out


def merge(left_df: pd.DataFrame, right_df: pd.DataFrame, on: str, how: str = "outer",
          update: bool = True, deduplicate: bool = False, **kwarg):
    assert all(isinstance(i, pd.DataFrame) for i in [left_df, right_df])
    _left_df = left_df.copy()
    _right_df = right_df.copy()
    if update:
        _left_df.update(_right_df)
    df1_unique_columns = [i for i in _right_df.columns if i not in _left_df.columns]
    out = _left_df.merge(_right_df.loc[:, [on] + df1_unique_columns], on=on, how=how, **kwarg)
    if deduplicate:
        out = deduplicate_df_by_row_merging(out, on)
    return out


def count_column_sizes(df: pd.DataFrame):
    return df.apply(lambda y: max(y.map(lambda x: len(str(x)))))


def excel_to_dfs_dict(file: str, **kwargs):
    sheets = pd.ExcelFile(file).sheet_names
    out = dict()
    for sheet in sheets:
        try:
            out[sheet] = pd.read_excel(file, sheet_name=sheet, **kwargs)
        except Exception:
            traceback.print_exc()
            print(f"Cannot parse sheet '{sheet}' from file '{file}'")
    return out


def dfs_dict_to_excel(d: dict, file: str, **kwargs):
    """
    :param d: {sheet_name <str>: sheet_dataframe <pandas.DataFrame>}
    :param file: table.xlx
    :return:
    """
    import xlsxwriter
    w = pd.ExcelWriter(file, engine="xlsxwriter")
    for sheet, df in d.items():
        try:
            df.to_excel(w, index=False, sheet_name=sheet, **kwargs)
        except Exception:
            traceback.print_exc()
            print(f"Cannot save sheet '{sheet}' to file '{file}'")
    w.save()


def remove_longest_columns(df: pd.DataFrame, size: int = 32767):  # M$ Excel cell size limit
    max_lengths = count_column_sizes(df)
    return df.loc[:, max_lengths.loc[max_lengths < size].index]


def dwell_df_on_column(df: pd.DataFrame, column_name: str):
    if column_name not in df.columns:
        return df
    return df.loc[:, [column_name] + [i for i in df.columns if i != column_name]]
