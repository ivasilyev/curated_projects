#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
from meta.scripts.Utilities import Utilities


class SampleDataLine:
    is_valid = False

    def __init__(self, sample_name: str, sample_read_files: list):
        self.name = sample_name.strip()
        self.reads = sorted(Utilities.remove_empty_values(sample_read_files))
        self._validate_reads()

    def __eq__(self, other):
        return self.name == other.name

    def __lt__(self, other):
        return self.name < other.name

    def _validate_reads(self):
        c = 0
        for read_file in self.reads:
            if not Utilities.is_file_valid(read_file, False):
                print("Not found the raw read file: '{}'".format(read_file))
                c += 1
        self.is_valid = c == 0

    @staticmethod
    def parse(d: dict):
        return SampleDataLine(d["sample_name"], d["raw_reads"])

    def export(self):
        return dict(sample_name=self.name, sample_reads=self.reads)


class SampleDataArray:
    def __init__(self):
        self.lines = dict()

    def __len__(self):
        return len(self.lines.keys())

    def validate(self):
        o = []
        for line in self.lines:
            if line.is_valid():
                o.append(line)
            else:
                print("Some read files are missing!")
        self.lines = sorted(o)

    @staticmethod
    def parse(d: dict):
        arr = SampleDataArray()
        arr.lines = {k: SampleDataLine.parse(d[k]) for k in d.keys()}
        return arr

    @staticmethod
    def load(file: str):
        return SampleDataArray.parse(json.loads(Utilities.load_string(file)))

    @staticmethod
    def generate(pair_2d_array: list, regex: str = "(.+)_S[0-9]+_L[0-9]+_R[0-9]+_[0-9]+"):
        arr = SampleDataArray()
        for read_files in pair_2d_array:
            read_files = sorted(read_files)
            sample_name = Utilities.safe_findall(regex, os.path.basename(read_files[0]))
            arr.lines[sample_name] = SampleDataLine(sample_name, read_files)
        return arr

    def export(self):
        return {k: self.lines[k].export() for k in self.lines}

    def to_dataframe(self):
        import pandas as pd
        return pd.DataFrame(list(self.export().values()))

    def dump(self, file: str):
        Utilities.dump_string(json.dumps(self.export(), sort_keys=True, indent=4), file)
