#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from meta.scripts.Utilities import Utilities


class SampleDataLine:
    is_valid = False

    def __init__(self, sample_name: str, sample_read_files: list):
        self.name = sample_name.strip()
        self.reads = Utilities.remove_empty_values(sample_read_files)
        self._validate_reads()

    def __eq__(self, other):
        return self.name == other.name

    def __lt__(self, other):
        return self.name < other.name

    def _validate_reads(self):
        c = 0
        for read_file in self.reads:
            if not os.path.isfile(read_file):
                print("Not found the raw read file: '{}'".format(read_file))
                c += 1
        self.is_valid = c == 0

    def export(self):
        return "\t".join([self.name, ";".join(self.reads)])


class SampleDataArray:
    def __init__(self):
        self.lines = []

    def __len__(self):
        return len(self.lines)

    def validate(self):
        o = []
        for line in self.lines:
            if line.is_valid():
                o.append(line)
            else:
                print("Some raw read files are missing!")
        self.lines = sorted(o)

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

    def to_dataframe(self):
        from io import StringIO
        return Utilities.load_tsv(StringIO(self.export()))
