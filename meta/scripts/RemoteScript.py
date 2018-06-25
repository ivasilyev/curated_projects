#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess


class RemoteScript:
    def __init__(self, url, output_dir):
        self._url = url
        if not output_dir:
            output_dir = os.getcwd()
        output_dir = (output_dir + "/", output_dir)[output_dir.endswith("/")]
        script_name = url.split("/")[-1]
        self.file = output_dir + script_name
        subprocess.getoutput("curl -fsSL {a} -o {b}".format(a=self._url, b=self.file))
