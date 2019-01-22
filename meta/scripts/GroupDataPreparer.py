#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd


class GroupDataPreparer:
    def __init__(self, raw_groupdata_file, processed_prefix, processed_suffix):
        self._df = pd.read_table(filepath_or_buffer=raw_groupdata_file,
                                 sep='\t',
                                 header='infer',
                                 names=["sample_name", "group_name"],
                                 engine='python')
        self._df["file_name"] = self._df.loc[:, "sample_name"].map(lambda x: processed_prefix + x + processed_suffix)
        self._df["file_exists"] = self._df.loc[:, "file_name"].map(lambda x: os.path.isfile(x))
        self._df = self._df.loc[self._df.loc[:, "file_exists"] == True]
        self.groups_string = "_".join(sorted(list(set(self._df.loc[:, "group_name"].values.tolist()))))
        if sum(self._df.loc[:, "file_exists"].values.tolist()) == 0:
            raise ValueError("Cannot find files by the mask: '{a}<sample name>{b}'".format(a=processed_prefix, b=processed_suffix))
    def get_groupdata(self):
        return self._df.loc[:, ["file_name", "group_name"]]
    def export_groupdata(self, output_file):
        self.get_groupdata().to_csv(path_or_buf=output_file, sep='\t', header=False, index=False)


class GroupDataAssemblyGuideLiner:
    def __init__(self, groupdata, prefix, suffix, index_column, value_column, output_dir):
        self.internal_launch_command = """python3 /home/docker/scripts/statistical_tools/groupdata2statistics.py \\
-g {GROUPDATA} \\
-p {PREFIX} \\
-s {SUFFIX} \\
-i {INDEX} \\
-v {VALUE} \\
-o {OUTPUT_DIR}""".format(GROUPDATA=groupdata, PREFIX=prefix, SUFFIX=suffix, INDEX=index_column, VALUE=value_column,
                          OUTPUT_DIR=output_dir)
        self.external_launch_command = """# Pre-setup to launch from different node for group data file \"{GROUPDATA}\" and value column \"{VALUE}\":

export IMG=ivasilyev/curated_projects:latest && \\
docker pull $IMG && \\
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it $IMG {CMD}

""".format(a=groupdata, GROUPDATA=groupdata, VALUE=value_column, CMD=self.internal_launch_command)
