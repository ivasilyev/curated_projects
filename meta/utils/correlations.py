#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import warnings
import numpy as np
import pandas as pd
from scipy import stats
from time import perf_counter
from itertools import combinations
from joblib import delayed, Parallel
from meta.utils.date_time import count_elapsed_seconds


def _process_out(__d: dict):
    _d = dict(__d)
    _d["significance_level"] = sum([_d["p_value"] < i for i in (0.01, 0.05, 0.1)])
    _d["denoted_correlation"] = "{}{}".format(
        round(_d["correlation_value"], 2),
        "*" * _d["significance_level"]
    )
    # Based on the Chaddock's correlation scale
    _d["chaddock_tightness"] = sum([
        _d["correlation_value"] >= i
        for i in (0.1, 0.3, 0.5, 0.7, 0.9)
    ])
    return _d


def mp_correlation_count(sub_df: pd.DataFrame, correlation_function = stats.spearmanr):
    """
    :param sub_df: pandas.DataFrame
    :param correlation_function: callable
    :return: dict
    """
    d = dict()
    if sub_df.shape[1] > 2:
        raise ValueError(f"Too much columns: {sub_df.columns}")
    d["x"], d["y"] = sub_df.columns
    d.update(dict(
        is_correlation_valid=False,
        p_value=0,
        correlation_value=0,
    ))
    d.update(_process_out(d))
    if len(set(sub_df.columns)) < 2:
        return d
    d["correlation_value"], d["p_value"] = correlation_function(
        sub_df.iloc[:, 0].values,
        sub_df.iloc[:, 1].values
    )
    if np.isnan(d["correlation_value"]) or sub_df.sum().sum() == 0:
        return d
    d["is_correlation_valid"] = True
    d.update(_process_out(d))
    return d


def slice_and_correlate_df(df: pd.DataFrame, correlation_function = stats.spearmanr):
    queue = [
        dict(
            sub_df=df.loc[:, i],
            correlation_function=correlation_function
        ) for i in combinations(set(df.columns), 2)
    ]
    now = perf_counter()
    print("Correlations to count: {}.".format(len(queue)))
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        correlations = Parallel(n_jobs=-1)(
            delayed(mp_correlation_count)(**i) for i in queue
        )
    print(f"Correlation count completed in {count_elapsed_seconds(now)}")
    correlation_total_df = pd.DataFrame(correlations)
    correlation_significant_df = correlation_total_df.query(
        "significance_level > 0 and "
        "chaddock_tightness > 0 and "
        "is_correlation_valid == True"
    ).sort_values("correlation_value", ascending=False)
    return dict(
        correlations_total=correlation_total_df,
        correlations_significant=correlation_significant_df,
    )
