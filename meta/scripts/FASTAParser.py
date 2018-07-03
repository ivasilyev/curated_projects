#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from meta.scripts.utilities import remove_empty_values


class FASTA:
    """
    This class is an attempt to apply NCBI standards to single FASTA.
    Consumes one header followed by sequence.
    """
    def __init__(self, single_fasta):
        self._body = re.sub("\n+", "\n", single_fasta.replace('\r', ''))
        try:
            self.header = re.findall("^>(.+)", self._body)[0].strip()
            self.sequence = "\n".join(self.chunk_string(re.sub("[^A-Za-z]", "", self._body.replace(self.header, "")), 70)).strip().upper()
        except IndexError:
            raise ValueError("Cannot parse the header for sequence: {}".format(self._body))
        # Nucleotide sequence has only AT(U)GC letters. However, it may be also protein FASTA.
    @staticmethod
    def chunk_string(string, length):
        return [string[0 + i:length + i] for i in range(0, len(string), length)]
    def __len__(self):
        return len(self.sequence)
    def to_dict(self):
        return {self.header: self.sequence}
    def to_str(self):
        return "\n".join([">" + self.header, self.sequence])
    def set_header(self, header):
        if header:
            if len(header) > 0:
                self.header = header


class FASTAParser:
    """
    This class parses FASTA file read as single string
    """
    def __init__(self, fastas_string):
        self._fastas_string = fastas_string
        self._raw_fastas_list = [">{}".format(j) if not j.startswith(">") else j for j in [i.strip() for i in re.split("\n>", self._fastas_string)]]
        self._parsed_fastas_list = remove_empty_values([FASTA(i) for i in self._raw_fastas_list])
    def get_full_length(self):
        return sum([len(i) for i in self._parsed_fastas_list])
    def get_full_sequence(self):
        return "\n".join([i.to_str() for i in self._parsed_fastas_list])
    def get_fastas_list(self):
        return self._parsed_fastas_list
    def get_headers_list(self):
        return [i.header for i in self._parsed_fastas_list]
