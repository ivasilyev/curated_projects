#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A worker script
"""

import os
import sys
import xlrd
import joblib
import warnings
import numpy as np
import pandas as pd
from time import sleep
from scipy import stats
from itertools import combinations
from meta.scripts.Utilities import Utilities
from meta.scripts.utils.pandas_utils import load_tsv, dump_tsv
from ashestopalov.nutrition.obesity_metagenomes.ProjectDescriber import ProjectDescriber


def mp_correlation_count(t: tuple):
    d = dict(columns=t, correlation=1, p_value=0, denoted_correlation="1", significance_level=0)
    if t[0] == t[-1]:
        return d
    d["correlation"], d["p_value"] = stats.spearmanr(correlation_df.loc[:, t])
    if np.isnan(d["correlation"]):
        return d
    d["significance_level"] = sum([d["p_value"] < i for i in (0.01, 0.05, 0.1)])
    d["denoted_correlation"] = "{}{}".format(round(d["correlation"], 2), "*" * d["significance_level"])
    return d


sleep(np.random.randint(90))

remote_queue = os.path.join(ProjectDescriber.DATA_DIR, "correlation_data", "group_datasets", "tables.txt")
correlation_tables = Utilities.remove_empty_values(Utilities.load_list(remote_queue))
if len(correlation_tables) == 0:
    print("Empty remote queue")
    sys.exit(0)

Utilities.dump_list(correlation_tables[1:], remote_queue)

correlation_table = correlation_tables[0]
group_name = os.path.splitext(os.path.basename(correlation_table))[0].replace("dataset_", "")
print("Now processing: '{}'".format(group_name))

correlation_df = load_tsv(correlation_table)
queue = list(combinations(correlation_df.columns, 2))

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    correlations = joblib.Parallel(n_jobs=-1)(joblib.delayed(mp_correlation_count)(i) for i in queue)

out_dir = os.path.join(ProjectDescriber.DATA_DIR, "correlation_data", "group_results", group_name)


Utilities.dump_list(correlations, os.path.join(out_dir, "all_results_for_{}.txt".format(group_name)))

valid_correlations = [dict(
    pair=" vs ".join(i["columns"]), correlation=i["correlation"], p_value=i["p_value"],
    denoted_correlation=i["denoted_correlation"], significance_level=i["significance_level"])
                      for i in correlations if not np.isnan(i["correlation"])]

valid_correlation_df = pd.DataFrame(valid_correlations)

for significance_level in sorted(set(valid_correlation_df["significance_level"].values)):
    significant_df = valid_correlation_df.query("significance_level == {}".format(significance_level))
    dump_tsv(significant_df, os.path.join(out_dir, "results_with_significance_{}_for_{}.tsv".format(
        significance_level, group_name)))
