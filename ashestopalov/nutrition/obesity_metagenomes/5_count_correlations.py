#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A worker script
"""

import os
import sys
import joblib
import warnings
import numpy as np
import pandas as pd
from time import sleep
from scipy import stats
from itertools import combinations, product
from meta.scripts.Utilities import Utilities
from meta.scripts.utils.pandas_utils import load_tsv, dump_tsv
from ashestopalov.nutrition.obesity_metagenomes.ProjectDescriber import ProjectDescriber


def mp_correlation_count(t: tuple):
    def _process_out():
        d["significance_level"] = sum([d["p_value"] < i for i in (0.01, 0.05, 0.1)])
        d["denoted_correlation"] = "{}{}".format(round(d["spearman_correlation"], 2),
                                                 "*" * d["significance_level"])
        # Based on the Chaddock's correlation scale
        d["chaddock_tightness"] = sum([d["spearman_correlation"] >= i
                                       for i in (0.1, 0.3, 0.5, 0.7, 0.9)])
        return d

    d = dict(feature_1=t[0], feature_2=t[-1], spearman_correlation=0, p_value=0,
             is_correlation_valid=False)
    d = _process_out()
    if t[0] == t[-1]:
        return d
    sub_df = correlation_df.loc[:, t]
    d["spearman_correlation"], d["p_value"] = stats.spearmanr(sub_df)
    if np.isnan(d["spearman_correlation"]) or sub_df.sum().sum() == 0:
        return d
    d["is_correlation_valid"] = True
    return _process_out()


try:
    print("Running on the node {}".format(os.uname()[1]))
except:
    pass

sleep(np.random.randint(90))
print("Polling the queue")

remote_queue = os.path.join(ProjectDescriber.DATA_DIR, "correlation_data", "group_datasets", "tables.txt")
correlation_tables = Utilities.remove_empty_values(Utilities.load_list(remote_queue))
if len(correlation_tables) == 0:
    print("Empty remote queue")
    sys.exit(0)

Utilities.dump_list(correlation_tables[1:], remote_queue)

correlation_table = correlation_tables[0]
print("Now processing: '{}'".format(correlation_table))

group_name = os.path.splitext(os.path.basename(correlation_table))[0]
out_dir = os.path.join(ProjectDescriber.DATA_DIR, "correlation_data", "group_results", group_name)

correlation_df = load_tsv(correlation_table).dropna(axis=0, how="any")
feature_groups = sorted(set([i.split("@")[0] for i in correlation_df.columns]))

if len(feature_groups) < 2:
    queue = list(combinations(correlation_df.columns, 2))
else:
    queue = list(product(*[[j for j in correlation_df.columns if j.startswith("{}@".format(i))]
                           for i in feature_groups]))

post_correlation_table = os.path.join(out_dir, "all_results_for_{}.tsv".format(group_name))
valid_correlation_table = os.path.join(out_dir, "somewhat_significant_results_for_{}.tsv".format(
    group_name))
if all(os.path.isfile(i) for i in [post_correlation_table, valid_correlation_table]):
    print("The output files exist: '{}', '{}'".format(post_correlation_table,
                                                      valid_correlation_table))
    existing_correlation_df = load_tsv(post_correlation_table)
    if existing_correlation_df.shape[0] == len(queue):
        print("Exiting, the existing table seems to be already fully processed: '{}'".format(
            post_correlation_table))
        sys.exit(0)

print("Starting the main processing, correlations to count: {}".format(len(queue)))
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    correlations = joblib.Parallel(n_jobs=-1)(joblib.delayed(mp_correlation_count)(i) for i in queue)


post_correlation_df = pd.DataFrame(correlations)
valid_correlation_df = post_correlation_df.query(
    "significance_level > 0 and chaddock_tightness > 0 and is_correlation_valid == True")

print("The valid table has {} rows".format(valid_correlation_df.shape[0]))

for df, table in zip([post_correlation_df, valid_correlation_df],
                     [post_correlation_table, valid_correlation_table]):
    dump_tsv(df.sort_values("spearman_correlation", ascending=False), table)

print("Processed '{}'".format(correlation_table))
