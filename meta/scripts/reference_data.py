#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
from abc import ABC, ABCMeta, abstractmethod
from meta.scripts.guidelines import dump_index_guide
from meta.utils.io import dump_dict, load_dict
from meta.utils.pandas import dump_tsv, load_tsv
from meta.utils.primitive import object_to_dict
from meta.utils.file_system import backup_file, find_file_by_tail, is_file_valid


class ReferenceData:
    def __init__(self):
        self.refdata_file = ""
        self.refdata_dict = dict()

    @property
    def sequences(self):
        return tuple(sorted([i for i in self.refdata_dict.keys() if i.startswith("sequence_")]))

    def parse(self, refdata_file):
        self.refdata_file = refdata_file
        self.refdata_dict = load_dict(self.refdata_file)
        print(f"Parsed {len(self.sequences)} reference sequences")

    def get_sequence_dict(self, number: int = 1):
        n = len(self.sequences)
        if n == 0:
            raise ValueError("This ReferenceData instance does not have any sequences")
        if number > n:
            raise ValueError(f"Too much number: '{number}' ({n} max)")
        return self.refdata_dict[self.sequences[number - 1]]

    @staticmethod
    def load(refdata_file: str):
        refdata = ReferenceData()
        refdata.parse(refdata_file)
        return refdata

    @staticmethod
    def find_and_load_refdata(directory: str):
        print(f"Looking for '*refdata.json' in '{directory}'")
        refdata_file = find_file_by_tail(directory, "_refdata.json")
        if len(refdata_file) == 0:
            print(f"Cannot find a ReferenceData compatible file within the directory '{directory}'")
            raise ValueError
        print(f"Reference data found at '{refdata_file}'")
        return ReferenceData.load(refdata_file)

    def update_metadata(self, d: dict):
        if "metadata" not in self.refdata_dict.keys():
            self.refdata_dict["metadata"] = dict()
        self.refdata_dict["metadata"].update(d)

    def dump(self, backup=False):
        if backup:
            _ = backup_file(self.refdata_file)
        dump_dict(self.refdata_dict, self.refdata_file)
        print(f"Reference Data saved to '{self.refdata_file}'")


class AnnotatorTemplate(ABC):
    __metaclass__ = ABCMeta
    INDEX_COL_NAME = "reference_id"

    def __init__(self):
        super().__init__()
        self.refdata = ReferenceData()
        self.annotation_file = ""
        self.annotation_df = pd.DataFrame()
        self._annotation_df = pd.DataFrame()

    def load(self):
        # Most of reference sequences are short enough to not be split
        self.annotation_file = self.refdata.get_sequence_dict()["annotation"]
        print(f"Using the annotation file: '{self.annotation_file}'")
        self.annotation_df = load_tsv(self.annotation_file)
        print(f"Loaded annotation file with shape '{self.annotation_df.shape}'")
        self._annotation_df = self.annotation_df.copy()

    def load_refdata(self, refdata_file: str):
        self.refdata.load(refdata_file)

    def validate(self):
        values = sorted(self.annotation_df[self.INDEX_COL_NAME].values.tolist())
        if len(values) != len(set(values)):
            raise ValueError(f"The index '{self.INDEX_COL_NAME}' is not uniquely valued!")
        _values = sorted(self._annotation_df[self.INDEX_COL_NAME].values.tolist())
        if values != _values:
            raise ValueError(f"Excessive values: {np.setxor1d([values, _values])}")

    def dump(self):
        backup = backup_file(self.annotation_file)
        self.refdata.dump()
        print(f"Annotation file backup saved to: '{backup}'")
        dump_tsv(df=self.annotation_df.sort_values(self.INDEX_COL_NAME),
                 table_file=self.annotation_file)
        print(f"Annotation file saved to: '{backup}'")

    def annotate(self):
        return

    def run(self):
        self.load()
        self.annotate()
        self.dump()


class ReferenceDescriberTemplate(ABC):
    __metaclass__ = ABCMeta
    NAME = ""
    VERSION = ""
    DESCRIPTION = ""
    DOCUMENTATION = ""
    WEBSITE = ""
    REFDATA_FILE = ""

    def __init__(self):
        super().__init__()
        self._refdata_parser = ""

    @property
    def ALIAS(self):
        return "{}_v{}".format(*[i.lower() for i in (self.NAME, self.VERSION)])

    def to_dict(self):
        return object_to_dict(self)

    def describe(self):
        return

    def parse_args(self):
        from argparse import ArgumentParser, RawTextHelpFormatter
        parser = ArgumentParser(formatter_class=RawTextHelpFormatter,
                                description=f"Describe and annotate {self.NAME}",
                                epilog="")
        parser.add_argument("-o", "--output", metavar="<dir>", required=True,
                            help="Output directory")
        namespace = parser.parse_args()
        return namespace.output


class SequenceRetrieverTemplate(ABC):
    __metaclass__ = ABCMeta
    NUCLEOTIDE_FASTA = ""
    PROTEIN_FASTA = ""
    REFERENCE_ROOT_DIRECTORY = ""
    REFERENCE_ANNOTATION = ""
    REFDATA_FILE = ""

    def __init__(self, describer: ReferenceDescriberTemplate):
        super().__init__()
        self.refdata = ReferenceData()

        self._reference_describer = describer
        self._reference_describer.REFDATA_FILE = self.REFDATA_FILE
        self._reference_describer.VERSION = self.VERSION

    @property
    def REFERENCE_DOWNLOAD_DIRECTORY(self):
        return os.path.join(self.REFERENCE_ROOT_DIRECTORY, self._reference_describer.ALIAS)

    @property
    def REFERENCE_INDEX_DIRECTORY(self):
        return os.path.join(self.REFERENCE_DOWNLOAD_DIRECTORY, "index")

    @property
    def VERSION(self):
        return self._reference_describer.VERSION

    @VERSION.setter
    def VERSION(self, value: str):
        self._reference_describer.VERSION = value
        print(f"The latest {self._reference_describer.NAME} reference version is {value}")

    @VERSION.deleter
    def VERSION(self):
        self._reference_describer.VERSION = ""

    def reset_nucleotide_fasta(self):
        # `cook_the_reference.py` takes alias directly from the file name
        self.NUCLEOTIDE_FASTA = os.path.join(self.REFERENCE_DOWNLOAD_DIRECTORY,
                                             f"{self._reference_describer.ALIAS}.fna")

    def reset_protein_fasta(self):
        self.PROTEIN_FASTA = f"{os.path.splitext(self.NUCLEOTIDE_FASTA)[0]}.faa"

    def get_latest_version(self):
        return

    def download(self):
        return

    def retrieve(self):
        return

    def create_nucleotide_fasta_symlink(self, real_file: str, default_nfasta: bool = False):
        if default_nfasta:
            self.reset_nucleotide_fasta()
        os.makedirs(os.path.dirname(self.NUCLEOTIDE_FASTA), exist_ok=True)
        print(f"Create symlink: '{real_file}' <===> '{self.NUCLEOTIDE_FASTA}'")
        os.symlink(real_file, self.NUCLEOTIDE_FASTA)

    def pick_refdata(self):
        try:
            self.refdata = ReferenceData.find_and_load_refdata(self.REFERENCE_INDEX_DIRECTORY)
            return True
        except ValueError:
            if is_file_valid(self.NUCLEOTIDE_FASTA):
                dump_index_guide(self.NUCLEOTIDE_FASTA, self.REFERENCE_INDEX_DIRECTORY)
            return False
