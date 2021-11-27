#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from abc import ABC, ABCMeta, abstractmethod
from meta.scripts.Utilities import Utilities


class ReferenceData:
    def __init__(self):
        self.refdata_file = ""
        self.refdata_dict = dict()
        self.sequences = ()

    def load(self, refdata_file):
        self.refdata_file = refdata_file
        self.refdata_dict = Utilities.load_dict(self.refdata_file)
        self.sequences = tuple(sorted([i for i in self.refdata_dict.keys if i.startswith("sequence_")]))

    def get_sequence_dict(self, number: int = 1):
        n = len(self.sequences)
        if number > n:
            raise ValueError(f"Too much number: '{number}' ({n} max)")
        return self.refdata_dict[self.sequences[number - 1]]

    @staticmethod
    def create(refdata_file: str):
        refdata = ReferenceData()
        refdata.load(refdata_file)
        return refdata

    @staticmethod
    def find_refdata(directory: str):
        refdata_file = Utilities.find_file_by_tail(directory, "_refdata.json")
        if len(refdata_file) == 0:
            raise ValueError(f"Cannot find a RefData file within the directory '{directory}'")
        return ReferenceData.create(refdata_file)

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

    def load(self, refdata_file: str):
        self.refdata.load(refdata_file)
        # Most of reference sequences are short enough to not be split
        self.annotation_file = self.refdata.get_sequence_dict()["annotation"]
        self.annotation_df = Utilities.load_tsv(self.annotation_file)

    def dump(self):
        _ = Utilities.backup_file(self.annotation_file)
        self.refdata.dump()
        Utilities.dump_tsv(self.annotation_df, self.annotation_file)

    def annotate(self):
        return


class SequenceRetrieverTemplate(ABC):
    __metaclass__ = ABCMeta

    def __init__(self):
        super().__init__()

    def download(self):
        return

    def retrieve(self):
        return


class ReferenceDescriberTemplate(ABC):
    __metaclass__ = ABCMeta
    NAME = ""
    VERSION = ""
    DESCRIPTION = ""
    DOCUMENTATION = ""
    WEBSITE = ""
    REFDATA = ""

    def __init__(self):
        super().__init__()
        self._refdata_parser = ""

    @property
    def ALIAS(self):
        return "{}_v{}".format(*[i.lower() for i in (self.NAME, self.VERSION)])

    def to_dict(self):
        return Utilities.object_to_dict(self)

    @staticmethod
    def get_index_guide(raw_nfasta_file):
        from meta.scripts.LaunchGuideLiner import LaunchGuideLiner
        index_directory = os.path.join(os.path.dirname(raw_nfasta_file), "index")
        LaunchGuideLiner.get_index_guide(
            index_directory=index_directory,
            raw_nfasta_file=raw_nfasta_file)
        return index_directory

    def describe(self):
        return
