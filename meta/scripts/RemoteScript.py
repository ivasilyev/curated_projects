#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess


class RemoteScript:
    def __init__(self, url, output_dir=os.getcwd()):
        self._url = url
        output_dir = (output_dir + "/", output_dir)[output_dir.endswith("/")]
        script_name = url.split("/")[-1]
        self.file = output_dir + script_name
    def download(self):
        subprocess.getoutput("curl -fsSL {a} -o {b}".format(a=self._url, b=self.file))
