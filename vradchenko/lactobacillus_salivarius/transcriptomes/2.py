#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%

import os
import pandas as pd
from meta.scripts.Utilities import Utilities

#%%

sra_dir = "/data1/bio/projects/vradchenko/lactobacillus_salivarius/sra"
sra_df = Utilities.load_tsv(os.path.join(sra_dir, "sra.tsv"))

queue = [{"func": Utilities.count_reads_statistics,
          "kwargs": {"reads_file": i, "type_": "fastq_gz"}}
         for i in Utilities.scan_whole_dir(os.path.join(sra_dir, "reads"))]

raw_reads_base_stats = Utilities.multi_core_queue(Utilities.wrapper, queue, async_=True)

#%%

raw_reads_base_stat_df = pd.DataFrame(raw_reads_base_stats)
raw_reads_base_stat_df["reads_file"] = raw_reads_base_stat_df["reads_file"].apply(os.path.basename)
raw_reads_base_stat_df["sample_name"] = raw_reads_base_stat_df["reads_file"].str.extract(r"(.+)\[")

Utilities.dump_tsv(raw_reads_base_stat_df, os.path.join(sra_dir, "raw_reads_base_stats.tsv"))
