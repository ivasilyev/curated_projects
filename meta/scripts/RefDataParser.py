#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class RefDataLine:
    """
    The class describes a single REFDATA dictionary
    """
    def __init__(self, parsed_dictionary):
        self.bowtie_index_mask = parsed_dictionary["ebwt_mask"]
        self.bowtie2_index_mask = parsed_dictionary["bt2_mask"]
        self.samtools_index_file = parsed_dictionary["fai"]
        self.bedtools_genome_file = parsed_dictionary["genome"]
        self.annotation_file = parsed_dictionary["annotation"]


class RefDataParser:
    """
    The class analyzes JSONs or REFDATA, tab-delimited text files containing paths to reference index files.
    REFDATA format:
    <Reference DNA fasta chunk 1> <Bowtie index mask> <Bowtie2 index mask> <SamTools index> <BedTools genome length> <Annotation file>
    <Reference DNA fasta chunk 2> <Bowtie index mask> <Bowtie2 index mask> <SamTools index> <BedTools genome length> <Annotation file>
    ...
    JSON keys (corresponding to above):
    {'sequence_1': {'reference_nfasta': '', 'ebwt_mask': '', 'bt2_mask': '', 'fai': '', 'genome': '', 'annotation': ''},
    'sequence_2': {'reference_nfasta': '', 'ebwt_mask': '', 'bt2_mask': '', 'fai': '', 'genome': '', 'annotation': ''},
    ...}
    """
    def __init__(self, reference_data_file_name):
        self._file_wrapper = open(reference_data_file_name, mode="r", encoding="utf-8")
        self._refdata_keys_list = self.get_refdata_keys_list()
        if reference_data_file_name.endswith(".json"):
            self._body_dict = self._parse_json_refdata()
        else:
            self._body_dict = self._parse_table_refdata()
        try:
            self._verify_json_refdata(self._body_dict)
        except ValueError:
            print("Bad reference data file: {}".format(reference_data_file_name))
            raise
        self.refdata_lines_dict = {k: RefDataLine(self._body_dict[k]) for k in self._body_dict}

    @staticmethod
    def get_refdata_keys_list():
        return ["reference_nfasta", "ebwt_mask", "bt2_mask", "fai", "genome", "annotation"]

    def _parse_json_refdata(self):
        import json
        return json.load(self._file_wrapper)

    def _parse_table_refdata(self):
        d = {}
        counter = 1
        for line in self._file_wrapper:
            if len(line.strip()):
                d["sequence_{}".format(counter)] = {k: v for k, v in zip(self._refdata_keys_list, [i.strip() for i in line.split("\t")])}
        return d

    def _verify_json_refdata(self, d: dict):
        import os
        if len(d) == 0:
            raise ValueError("Empty sample data!")
        for k1 in d:
            if len([i for i in list(d) if i == k1]) > 1:
                raise ValueError("Repeating key: {}".format(k1))
            for k_r in [i for i in self._refdata_keys_list if i != "reference_nfasta"]:
                if k_r not in list(d[k1]):
                    print("Warning! Missing value for the key: '{}'".format(k_r))
            for k2 in d[k1]:
                f = d[k1][k2]
                if k2 not in ["alias", "reference_nfasta", "ebwt_mask", "bt2_mask"] and not os.path.isfile(f):
                    print("Warning! Not found file: '{}', keys: '{}', '{}'".format(f, k1, k2))

    def get_parsed_list(self):
        return [self.refdata_lines_dict[k] for k in self.refdata_lines_dict]

    def get_refdata_line_by_index(self, idx: int):
        return self.get_parsed_list()[idx]

    def get_refdata_line_by_key(self, key: str):
        return self.refdata_lines_dict.get(key)
