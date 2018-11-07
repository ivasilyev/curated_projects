#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class ReferenceDescriberTemplate():
    def __init__(self, alias: str, name: str, description: str, documentation: str):
        self.alias = alias
        self.name = name
        self.description = description
        self.documentation = documentation
    def export(self):
        print("""
class ReferenceDescriber:
    alias = "{ALIAS}"
    name = "{NAME}"
    description = "{DESCRIPTION}"
    documentation = "{DOCUMENTATION}"
    refdata = "PASTE_HERE_AFTER_INDEXING"
""".format(ALIAS=self.alias, NAME=self.name, DESCRIPTION=self.description, DOCUMENTATION=self.documentation))
    @staticmethod
    def parse_refdata(refdata):
        from meta.scripts.RefDataParser import RefDataParser
        refdata_parser = RefDataParser(refdata)
        return refdata_parser.get_parsed_list()
