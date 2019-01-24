#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from meta.scripts.Utilities import Utilities


class PivotSplitter:
    def __init__(self, pivot_df: pd.DataFrame, value_col_name: str, ):
        self.pivot_df = pivot_df
        self.value_col_name = value_col_name
        self._sample_names_list = []
    def split(self, output_dir: str):
        output_dir = Utilities.ends_with_slash(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        # Note: the dataframe must have only index and value columns
        for sample_col_name in list(self.pivot_df):
            sample_name = Utilities.filename_only(sample_col_name).split("_")[0]
            sample_file_name = "{}{}.tsv".format(output_dir, sample_name)
            self.pivot_df[sample_col_name].reset_index().rename(columns={sample_col_name: self.value_col_name}).to_csv(
                sample_file_name, sep="\t", header=True, index=False)
            self._sample_names_list.append(sample_file_name)
    def get_groupdata(self, group_name: str):
        return pd.DataFrame([{"sample_name": i, "group_name": group_name} for i in self._sample_names_list]).loc[:, ["sample_name", "group_name"]]
