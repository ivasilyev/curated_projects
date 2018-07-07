#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import yaml
import jinja2


class ChartGenerator:
    def __init__(self, config, template):
        self._raw_config = self._read_file_or_url(config)
        self._raw_template = self._read_file_or_url(template)
    @staticmethod
    def _read_file_or_url(path):
        if os.path.isfile(path):
            open(file=path, mode="r", encoding="utf-8")
            with open(path, 'r') as f:
                return f.read()
        else:
            if path.startswith("http") or path.startswith("www."):
                return subprocess.getoutput("curl -fsSL {}".format(path))
            else:
                raise ValueError("Cannot load  file or URL: '{}'".format(path))
    def _template2yaml(self, buffered_template):
        t = jinja2.Template(buffered_template)
        return t.render(yaml.load(self._raw_config))
    def export_yaml(self, file):
        s = self._template2yaml(self._raw_template)
        with open(file=file, mode="w", encoding="utf-8") as f:
            f.write(s)
