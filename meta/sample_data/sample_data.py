#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
import os
import re
import json
from meta.scripts.Utilities import Utilities
from argparse import ArgumentParser, RawTextHelpFormatter


class SampleDataLine:
    def __init__(self, sample_name: str, sample_read_files: list):
        self.state = dict()
        self.name = sample_name.strip()
        self.reads = sorted(Utilities.remove_empty_values(sample_read_files))
        self.taxa = ""
        self.is_valid = False
        self._validate_reads()

    def __eq__(self, other):
        return self.name == other.name

    def __lt__(self, other):
        return self.name < other.name

    def __len__(self):
        return len(self.reads)

    def __repr__(self):
        return "SampleDataLine with name '{}' and {} reads".format(self.name, len(self))

    def _validate_reads(self):
        c = 0
        for read_file in self.reads:
            if not Utilities.is_file_valid(read_file, False):
                print("Not found the raw read file: '{}'".format(read_file))
                c += 1
        self.is_valid = c == 0

    @staticmethod
    def parse(d: dict):
        return SampleDataLine(d["name"], d["reads"])

    def export(self):
        d = dict(name=self.name, reads=self.reads, taxa=self.taxa)
        d.update(self.state)
        return d


class SampleDataArray:
    def __init__(self):
        self.lines = dict()

    def __len__(self):
        return len(self.lines.keys())

    def __repr__(self):
        return "SampleDataArray with {} lines".format(len(self))

    @property
    def values(self):
        return self.lines.values()

    def validate(self):
        d = dict()
        for key in self.lines:
            line = self.lines[key]
            if line.is_valid:
                d[key] = line
            else:
                print("Invalid sample data for the sample: '{}'".format(line.name))
        self.lines = d

    def update_lines_state(self, d: dict):
        _ = [i.state.update(d) for i in self.lines.values()]

    @staticmethod
    def parse(d: dict):
        arr = SampleDataArray()
        arr.lines = {k: SampleDataLine.parse(d[k]) for k in d.keys()}
        arr.validate()
        return arr

    @staticmethod
    def load(file: str):
        return SampleDataArray.parse(json.loads(Utilities.load_string(file)))

    @staticmethod
    def generate(pair_2d_array: list, regex: str = DEFAULT_REGEX,
                 extension: str = DEFAULT_READS_EXTENSION):
        arr = SampleDataArray()
        for sample_read_files in pair_2d_array:
            sample_read_files = sorted(sample_read_files)
            sample_file = os.path.basename(sample_read_files[0])
            sample_name = Utilities.safe_findall(regex, re.sub(f"{extension}$", "", sample_file))
            if len(sample_name) == 0:
                raise ValueError(f"Cannot process the file '{sample_file}' with the regex '{regex}'")
            if any(sample_name not in i for i in sample_read_files):
                raise ValueError(f"Some files from the list '{sample_read_files}' do not contain {sample_name} parsed by the regex '{regex}'")
            if sample_name in arr.lines.keys():
                print(f"Duplicate sample data line key, the regex check is considered: '{sample_name}'")
                c = 0
                sample_name_ = str(sample_name)
                while sample_name in arr.lines.keys():
                    c += 1
                    sample_name = "{}.{}".format(sample_name_, c)
            arr.lines[sample_name] = SampleDataLine(sample_name, sample_read_files)
        return arr

    @staticmethod
    def generate_from_directory(directory: str, regex: str = DEFAULT_REGEX,
                                reads_extension: str = DEFAULT_READS_EXTENSION):
        pair_2d_array = Utilities.get_most_similar_word_pairs(
            Utilities.find_file_by_tail(directory, reads_extension))
        return SampleDataArray.generate(pair_2d_array, regex=regex, extension=reads_extension)

    def export(self):
        return {k: self.lines[k].export() for k in self.lines}

    def to_dataframe(self):
        import pandas as pd
        return pd.DataFrame(list(self.export().values()))

    def dump(self, file: str):
        Utilities.dump_dict(self.export(), file)
"""

"""
def parse_args():
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description="Generate sampledata based on directory".strip(),
        epilog=""
    )
    parser.add_argument("-i", "--input", nargs="+",
                        help="Input directory (directories)")
    parser.add_argument("-e", "--extension", default=DEFAULT_READS_EXTENSION,
                        help="Extension of reads files")
    parser.add_argument("-r", "--regex", default=DEFAULT_REGEX,
                        help="Regular expression to extract sample name(s) from file name(s) without extension")
    parser.add_argument("-t", "--taxa", default="",
                        help="(Optional) Taxonomy to add to all samples")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file")
    _namespace = parser.parse_args()
    return _namespace.input, _namespace.extension, _namespace.regex,  _namespace.taxa, _namespace.output


if __name__ == '__main__':
    inputDirs, inputExtension, inputRegex, inputTaxa, outputFile = parse_args()

    pair2dArray = []
    for input_dir in inputDirs:
        pair2dArray.extend(Utilities.get_most_similar_word_pairs(Utilities.find_file_by_tail(
            input_dir, inputExtension, multiple=True)))

    sampleDataArray = SampleDataArray.generate(pair2dArray, inputRegex)

    for sampleDataLine in sampleDataArray.values:
        sampleDataLine.taxa = inputTaxa

    sampleDataArray.dump(outputFile)
"""

from meta.utils.io import dump_dict
from meta.utils.language import tokenize_reads_file_name
from meta.utils.file_system import find_file_by_tail


DEFAULT_REGEX = "(.+).+S[0-9]+.*R[12]"
DEFAULT_READS_EXTENSION = "fastq.gz"


def create_sampledata_dict_from_list(reads_files: list):
    tokenized_reads_files = [
        tokenize_reads_file_name(i) for i in sorted(sorted(reads_files), key=len, reverse=True)
    ]
    out = dict()
    for token_dict in tokenized_reads_files:
        sample_name = token_dict["sample_name"]
        if sample_name in out.keys():
            out[sample_name]["reads"].append(token_dict["reads_file"])
        else:
            out[sample_name] = {
                "name": sample_name,
                "reads": [token_dict["reads_file"], ],
                "taxa": ""
            }
        out[sample_name]["reads"] = sorted(out[sample_name]["reads"])
    return out


def create_sampledata_dict_from_dir(directory: str, reads_extension: str = DEFAULT_READS_EXTENSION):
    reads_files = find_file_by_tail(
        directory, ".{}".format(reads_extension.strip(".")), multiple=True
    )
    return create_sampledata_dict_from_list(reads_files)


def parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description="Generate sample data based on scanned raw reads files in the directory".strip(),
        epilog=""
    )
    parser.add_argument("-i", "--input", nargs="+",
                        help="Input directory (directories)")
    parser.add_argument("-e", "--extension", default=DEFAULT_READS_EXTENSION,
                        help="Extension of reads files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file")
    _namespace = parser.parse_args()
    return _namespace.input, _namespace.extension, _namespace.output


if __name__ == '__main__':
    inputDirs, inputExtension, outputFile = parse_args()

    sampledata_dict = dict()
    for input_dir in inputDirs:
        sampledata_dict.update(create_sampledata_dict_from_dir(input_dir, inputExtension))

    dump_dict(sampledata_dict, outputFile)
    print(f"Sampledata successfully created in: '{outputFile}'")
