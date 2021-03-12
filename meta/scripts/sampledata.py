#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from meta.scripts.Utilities import Utilities


class SampleDataLine:
    def __init__(self, sample_name: str, sample_read_files: list):
        self.name = sample_name.strip()
        self.reads = Utilities.remove_empty_values(sample_read_files)

    def __eq__(self, other):
        return self.name == other.name

    def __lt__(self, other):
        return self.name < other.name

    def export(self):
        return "\t".join([self.name, ";".join(self.reads)])


class SampleDataArray:
    lines = []

    def __len__(self):
        return len(self.lines)

    def validate(self):
        self.lines = sorted([i for i in self.lines if i.exists])

    @staticmethod
    def generate(pair_2d_array: list, regex: str = "(.+)_S[0-9]+_L[0-9]+_R[0-9]+_[0-9]+"):
        arr = SampleDataArray()
        for raw_read_pair in pair_2d_array:
            raw_read_pair = sorted(raw_read_pair)
            basename = os.path.basename(raw_read_pair[0])
            sample_name = Utilities.safe_findall(regex, basename)
            arr.lines.append(SampleDataLine(sample_name, raw_read_pair))
        return arr

    def export(self):
        return "\n".join(["sample_name\traw_reads"] + [i.export() for i in self.lines])

    def to_tsv(self):
        from io import StringIO
        return Utilities.load_tsv(StringIO(self.export()))
