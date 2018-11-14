#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import pandas as pd
from collections import Counter
from meta.scripts.Utilities import Utilities


class CounterWrapper:
    @staticmethod
    def prepare_string(s: str):
        return re.sub(" +", " ", re.sub("[^a-z0-9\-]+", " ", s.lower())).strip()
    @staticmethod
    def count_words_in_series(series: pd.Series):
        lst = [i for i in CounterWrapper.prepare_string(" ".join(series.fillna("").values.tolist())).split(" ") if len(i) > 3]
        return Counter(lst)
    @staticmethod
    def dump_counter(counter: Counter, file: str):
        Utilities.dump_2d_array([("keyword", "occurrences")] + counter.most_common(), file=file)
