#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import ABC, ABCMeta, abstractmethod


class ReferenceDescriberTemplate(ABC):
    __metaclass__ = ABCMeta
    NAME = ""
    VERSION = ""
    ALIAS = ""
    DESCRIPTION = ""
    DOCUMENTATION = ""
    WEBSITE = ""
    REFDATA = ""

    def __init__(self):
        super().__init__()
        self._refdata_parser = ""

    def update_alias(self):
        self.ALIAS = "{}_v{}".format(self.NAME.lower(), self.VERSION.lower())
        return self.ALIAS

    def export(self):
        fields = """
    NAME = ""
    VERSION = ""
    ALIAS = ""
    DESCRIPTION = ""
    DOCUMENTATION = ""
    WEBSITE = ""
    REFDATA = ""
""".replace('""', '"{}"').format(self.NAME, self.VERSION, self.update_alias(), self.DESCRIPTION, self.DOCUMENTATION,
                                 self.WEBSITE, self.REFDATA)
        print("""Please update the following script lines: 

class ReferenceDescriber(ReferenceDescriberTemplate): 
{}""".format(fields))

    @staticmethod
    def get_index_guide(raw_nfasta_file):
        import os
        from meta.scripts.LaunchGuideLiner import LaunchGuideLiner
        LaunchGuideLiner.get_index_guide(
            index_directory=os.path.join(os.path.dirname(raw_nfasta_file), "index"),
            raw_nfasta_file=raw_nfasta_file)

    @staticmethod
    def find_refdata(index_dir):
        import os
        import subprocess
        cmd = 'ls -d {}/* | grep "_refdata.json"'.format(os.path.normpath(index_dir))
        out = subprocess.getoutput(cmd).strip()
        if out.count("\n") > 0:
            raise ValueError(
                "Cannot find single reference data file! \nPlease check out the shell command: `{}`".format(cmd))
        return out

    def get_refdata_dict(self):
        from meta.scripts.RefDataParser import RefDataParser
        self._refdata_parser = RefDataParser(self.REFDATA)
        return self._refdata_parser.refdata_lines_dict
