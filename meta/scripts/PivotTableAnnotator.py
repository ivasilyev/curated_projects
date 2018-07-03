#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd


class PivotTableAnnotator:
    def __init__(self, pivot_file, annotation_file):
        pivot_df = pd.read_table(pivot_file, sep='\t', header=0, engine='python')
        annotation_df = pd.read_table(annotation_file, sep='\t', header=0, engine='python')
        self.annotated_df = pd.merge(annotation_df, pivot_df, on="reference_id", how='outer')
    def export_annotated_pivot(self, output_file):
        self.annotated_df.to_csv(path_or_buf=output_file, sep='\t', header=True, index=False)
    def export_group_reports(self):
        pass

pivot_file = "/data1/bio/projects/dsafina/hp_checkpoints/TADB/groupdata2statistics/C_I_II_III_IV_srr_total_dataframe_annotated.tsv"
self_annotated_df = pd.read_table(pivot_file, sep='\t', header=0, engine='python')

import matplotlib.pyplot as plt
import os

multi_test = "fdr_bh"
single_test = "u-test"
collection_metric = "sum"
base_col_name = "reference_id"
annotation_col_name = "former_id"
output_dir = "/data1/bio/projects/dsafina/hp_checkpoints/TADB/groupdata2statistics/test/"
output_dir = (output_dir + "/", output_dir)[output_dir.endswith("/")]

bool_suffix = "_is_rejected_by_{m}_for_{s}".format(m=multi_test, s=single_test)
bool_col_names_list = [i for i in list(self_annotated_df) if i.endswith(bool_suffix)]
# for
bool_col_name = bool_col_names_list[0]
comparison_string = bool_col_name.replace(bool_suffix, "")
comparison_list = comparison_string.split("_vs_")
left_group_name = comparison_list[0]
right_group_name = comparison_list[1]
comparison_df = self_annotated_df.loc[self_annotated_df[bool_col_name] == True].set_index(base_col_name)
# if len(comparison_df) > 0:
left_group_pivot_df = comparison_df.loc[:, [annotation_col_name] + [i for i in comparison_df if i.split("_")[0] == left_group_name and i.split("_")[1].startswith("/")]]
right_group_pivot_df = comparison_df.loc[:, [annotation_col_name] + [i for i in comparison_df if i.split("_")[0] == right_group_name and i.split("_")[1].startswith("/")]]
comparison_dir = "{a}{b}/".format(a=output_dir, b=comparison_string)
os.makedirs(path=comparison_dir, exist_ok=True)
comparison_df.to_csv(path_or_buf="{}pivot.tsv".format(comparison_dir), sep='\t', header=True, index=True)
left_group_pivot_df.to_csv(path_or_buf="{a}raw_{b}.tsv".format(a=comparison_dir, b=left_group_name), sep='\t', header=True, index=True)
right_group_pivot_df.to_csv(path_or_buf="{a}raw_{b}.tsv".format(a=comparison_dir, b=right_group_name), sep='\t', header=True, index=True)
boxplots_dir = "{}boxplots/".format(comparison_dir)
# for base_id in comparison_df.index.tolist():
base_id = comparison_df.index.tolist()[0]
annotated_id = comparison_df.loc[base_id, annotation_col_name]






left_group_metric_df = comparison_df.loc[:, [annotation_col_name, "{g}_{m}".format(g=comparison_list[0], m=collection_metric)]]


list(comparison_df)



III_vs_IV_is_rejected_by_fdr_bh_for_h-test



