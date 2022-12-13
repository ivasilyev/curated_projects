#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from meta.utils.io import load_dict, dump_dict
from meta.sample_data.sample_data_line import SampleDataLine
from meta.utils.language import tokenize_reads_file_name
from meta.utils.file_system import find_file_by_tail


DEFAULT_READS_EXTENSION = "fastq.gz"


class SampleDataArray:
    def __init__(self, lines: dict = None):
        if lines is None:
            lines = dict()
        self.lines = lines
        self.is_valid = False
        self.validate()

    @property
    def names(self):
        return list(self.lines.keys())

    @property
    def values(self):
        return list(self.lines.values())

    def __len__(self):
        return len(self.names)

    def __repr__(self):
        return f"SampleDataArray with {len(self)} lines"

    @staticmethod
    def _is_line_valid(x: SampleDataLine):
        if x.is_valid:
            return True
        else:
            print(f"Invalid sample data for the sample: '{x.name}'")
            return False

    def validate(self):
        self.lines = {k: v for k, v in self.lines.items() if self._is_line_valid(v)}
        self.is_valid = len(self) > 0

    def update_lines_state(self, d: dict):
        _ = [i.state.update(d) for i in self.lines.values()]

    @staticmethod
    def import_from_dict(d: dict):
        arr = SampleDataArray({k: SampleDataLine.import_from_dict(v) for k, v in d.items()})
        arr.validate()
        return arr

    @staticmethod
    def import_from_dicts(x: list):
        out = dict()
        for d in x:
            line = SampleDataLine.import_from_dict(d)
            out[line.name] = line
        return SampleDataArray.import_from_dict(out)

    @staticmethod
    def load_dict(file: str):
        return SampleDataArray.import_from_dict(load_dict(file))

    @staticmethod
    def import_from_reads_files(reads_files: list):
        tokenized_reads_files = [
            tokenize_reads_file_name(i) for i in sorted(sorted(reads_files), key=len, reverse=True)
        ]
        arr = SampleDataArray()
        for token_dict in tokenized_reads_files:
            sample_name = token_dict["sample_name"]
            if sample_name in arr.names:
                line = arr.lines[sample_name]
                line.reads.append(token_dict["reads_file"])
            else:
                line = SampleDataLine(sample_name, [token_dict["reads_file"]])
            line.validate()
            arr.lines[sample_name] = line
        return arr

    @staticmethod
    def merge_arrays(arrays: list):
        out = SampleDataArray()
        for arr in arrays:
            out.lines.update(arr.lines)
        return out

    @staticmethod
    def import_from_dir(
        directory: str,
        reads_extension: str = DEFAULT_READS_EXTENSION
    ):
        reads_files = find_file_by_tail(
            directory, ".{}".format(reads_extension.strip(".")),
            multiple=True
        )
        return SampleDataArray.import_from_reads_files(reads_files)

    def export_to_dict(self):
        return {k: v.export_to_dict() for k, v in self.lines.items()}

    def to_dataframe(self):
        return pd.DataFrame(list(self.export_to_dict().values())).loc[:, ["sample_name", "reads"]]

    def dump_dict(self, file: str):
        self.validate()
        if self.is_valid:
            dump_dict(self.export_to_dict(), file)
            print(f"Exported sample data to: '{file}'")
        else:
            print("No sample data to export")

    def dump_table(self, file: str):
        self.validate()
        if self.is_valid:
            df = self.to_dataframe()
            os.makedirs(os.path.dirname(file), exist_ok=True)
            df.to_csv(file, encoding="utf-8", sep="\t", index=False, header=False)
            print(f"Exported sample data to: '{file}'")
        else:
            print("No sample data to export")
