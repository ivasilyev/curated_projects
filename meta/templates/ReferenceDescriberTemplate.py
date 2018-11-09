#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import ABC, ABCMeta, abstractmethod


class ReferenceDescriberTemplate(ABC):
    __metaclass__ = ABCMeta
    alias = ""
    name = ""
    description = ""
    documentation = ""
    refdata = ""
    def __init__(self):
        super().__init__()
    def export(self):
        print("""
class ReferenceDescriber(ReferenceDescriberTemplate):
    alias = "{ALIAS}"
    name = "{NAME}"
    description = "{DESCRIPTION}"
    documentation = "{DOCUMENTATION}"
    refdata = "{REFDATA}"
""".format(ALIAS=self.alias, NAME=self.name, DESCRIPTION=self.description, DOCUMENTATION=self.documentation,
           REFDATA=self.refdata))
    @staticmethod
    def parse_refdata(refdata):
        from meta.scripts.RefDataParser import RefDataParser
        refdata_parser = RefDataParser(refdata)
        return refdata_parser.get_parsed_list()
