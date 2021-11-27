#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from abc import ABC, ABCMeta, abstractmethod
from meta.scripts.Utilities import Utilities
from meta.scripts.guidelines import dump_index_guide


class ReferenceData:
    def __init__(self):
        self.refdata_file = ""
        self.refdata_dict = dict()
        self.sequences = ()

    def parse(self, refdata_file):
        self.refdata_file = refdata_file
        self.refdata_dict = Utilities.load_dict(self.refdata_file)
        self.sequences = tuple(sorted([i for i in self.refdata_dict.keys()
                                       if i.startswith("sequence_")]))

    def get_sequence_dict(self, number: int = 1):
        n = len(self.sequences)
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
        refdata_file = Utilities.find_file_by_tail(directory, "_refdata.json")
        if len(refdata_file) == 0:
            raise ValueError(f"Cannot find a RefData file within the directory '{directory}'")
        return ReferenceData.load(refdata_file)

    def update_metadata(self, d: dict):
        if "metadata" not in self.refdata_dict.keys():
            self.refdata_dict["metadata"] = dict()
        self.refdata_dict["metadata"].update(d)

    def dump(self, backup=False):
        if backup:
            _ = Utilities.backup_file(self.refdata_file)
        Utilities.dump_dict(self.refdata_dict, self.refdata_file)


class AnnotatorTemplate(ABC):
    __metaclass__ = ABCMeta

    def __init__(self):
        super().__init__()
        self.refdata = ReferenceData()
        self.annotation_file = ""
        self.annotation_df = pd.DataFrame()

    def parse_annotation(self):
        # Most of reference sequences are short enough to not be split
        self.annotation_file = self.refdata.get_sequence_dict()["annotation"]
        self.annotation_df = Utilities.load_tsv(self.annotation_file)

    def load_refdata(self, refdata_file: str):
        self.refdata.load(refdata_file)

    def dump(self):
        _ = Utilities.backup_file(self.annotation_file)
        self.refdata.dump()
        Utilities.dump_tsv(self.annotation_df, self.annotation_file)

    def annotate(self):
        return


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
        return Utilities.object_to_dict(self)

    def describe(self):
        return


class SequenceRetrieverTemplate(ABC):
    __metaclass__ = ABCMeta
    VERSION = ""
    NUCLEOTIDE_FASTA = ""
    REFERENCE_ROOT_DIRECTORY = ""
    REFERENCE_ANNOTATION = ""
    REFDATA_FILE = ""

    def __init__(self, describer: ReferenceDescriberTemplate):
        super().__init__()
        self._reference_describer = describer
        self._reference_describer.REFDATA_FILE = self.REFDATA_FILE
        self._reference_describer.VERSION = self.VERSION
        self.refdata = ReferenceData()

    @property
    def REFERENCE_FETCH_DIRECTORY(self):
        return os.path.join(self.REFERENCE_ROOT_DIRECTORY, self._reference_describer.ALIAS)

    @property
    def REFERENCE_INDEX_DIRECTORY(self):
        return os.path.join(self.REFERENCE_FETCH_DIRECTORY, "index")

    def download(self):
        return

    def retrieve(self):
        return

    def pick_refdata(self):
        try:
            self.refdata = ReferenceData.find_and_load_refdata(self.REFERENCE_INDEX_DIRECTORY)
            return True
        except ValueError:
            dump_index_guide(self.NUCLEOTIDE_FASTA, self.REFERENCE_INDEX_DIRECTORY)
            return False
