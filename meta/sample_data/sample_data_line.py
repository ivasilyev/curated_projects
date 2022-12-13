#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
from meta.utils.file_system import is_file_valid
from meta.utils.primitive import remove_empty_values


class SampleDataLine:
    def __init__(self, name: str, reads: list, taxa: dict = None):
        self.state = dict()
        self.name = name.strip()
        self.reads = remove_empty_values(reads)
        if taxa is None:
            taxa = dict()
        self.taxa = taxa
        self.is_valid = False
        self.validate()

    def __eq__(self, other):
        return self.name == other.name

    def __lt__(self, other):
        return self.name < other.name

    def __len__(self):
        return len(self.reads)

    def __repr__(self):
        return f"SampleDataLine with name '{self.name}' and {len(self)} reads"

    def _validate_name(self):
        self.is_valid = self.name is not None and len(self.name) > 0

    def _validate_reads(self):
        c = int(len(self.reads) == 0)
        out = list()
        for read_file in self.reads:
            if is_file_valid(file=read_file, report=False):
                out.append(read_file)
            else:
                print("Not found the reads file: '{}'".format(read_file))
                c += 1
        self.is_valid = c == 0
        self.reads = out

    def validate(self):
        self._validate_name()
        self._validate_reads()

    @staticmethod
    def import_from_dict(d: dict):
        """
        :param d: {
            'name': 'sample_name_1',
            'reads': ['reads_1', ...],
            'taxa':
                {
                    'genera': 'genera',
                    'species': 'species',
                    'strain': 'strain'
                }
            'state':
                {
                    'key': 'value'
                }
        }
        :return: SampleDataLine object
        """
        if any(i not in d.keys() for i in ["name", "reads"]):
            raise ValueError(f"Unable to parse: '{json.dumps(d)}'")
        out = SampleDataLine(d.get("name"), d.get("reads"))
        if "taxa" in d.keys():
            out.taxa = d.get("taxa")
        if "state" in d.keys():
            out.state = d.get("state")
        out.validate()
        return out

    def export_to_dict(self):
        d = dict(name=self.name, reads=self.reads, taxa=self.taxa, state=self.state)
        return d

