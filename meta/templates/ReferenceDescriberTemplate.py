#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import ABC, ABCMeta, abstractmethod


class ReferenceDescriberTemplate(ABC):
    __metaclass__ = ABCMeta
    alias = ""
    name = ""
    description = ""
    documentation = ""
    url = ""
    refdata = ""

    def __init__(self):
        super().__init__()
        self._refdata_parser = ""

    def export(self):
        print("""Please update the following script lines: 

class ReferenceDescriber(ReferenceDescriberTemplate):
    alias = "{ALIAS}"
    name = "{NAME}"
    description = "{DESCRIPTION}"
    documentation = "{DOCUMENTATION}"
    url = "{URL}"
    refdata = "{REFDATA}"

""".format(ALIAS=self.alias, NAME=self.name, DESCRIPTION=self.description, DOCUMENTATION=self.documentation,
           URL=self.url, REFDATA=self.refdata))

    def get_index_guide(self, raw_nfasta_file):
        import os
        from meta.scripts.LaunchGuideLiner import LaunchGuideLiner
        LaunchGuideLiner.get_index_guide(
            index_directory=os.path.join(os.path.dirname(raw_nfasta_file), "index", self.alias),
            raw_nfasta_file=raw_nfasta_file)

    def get_refdata_dict(self):
        from meta.scripts.RefDataParser import RefDataParser
        self._refdata_parser = RefDataParser(self.refdata)
        return self._refdata_parser.refdata_lines_dict
