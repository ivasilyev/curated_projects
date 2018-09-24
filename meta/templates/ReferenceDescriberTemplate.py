#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import ABC, ABCMeta, abstractmethod


class ReferenceDescriberTemplate(ABC):
    __metaclass__ = ABCMeta

    def __init__(self, alias: str, name: str, description: str, documentation: str, refdata: str):
        super().__init__()
        self.alias = alias
        self.name = name
        self.description = description
        self.documentation = documentation
        self.refdata = refdata

    @abstractmethod
    def export(self):
        print("""
    Database alias: {a}
    REFDATA linker: {b}
                  """.format(a=self.alias, b=self.refdata))

    @abstractmethod
    def parse_refdata(self):
        from meta.scripts.RefDataParser import RefDataParser
        return RefDataParser(self.refdata).get_parsed_list()
