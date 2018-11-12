#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${DOCKER_IMAGE_NAME} python3
"""

import os
import re
import subprocess
import pandas as pd
import numpy as np
from meta.templates.ReferenceDescriberTemplate import ReferenceDescriberTemplate


class ReferenceDescriber(ReferenceDescriberTemplate):
    alias = "vfdb_v2018.11.09"
    name = "VFDB"
    description = "A reference database for bacterial virulence factors"
    documentation = "http://www.mgc.ac.cn/VFs/main.htm"
    refdata = "/data/reference/VFDB/vfdb_v2018.11.09/index/vfdb_v2018.11.09_refdata.json"


class SequenceRetriever:
    # In the case of connection failure download data via chinese proxy from http://www.mgc.ac.cn/VFs/download.htm
    urls_dict = {"VFs description file": "http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz",
                 "Intra-genera VFs comparison tables": "http://www.mgc.ac.cn/VFs/Down/Comparative_tables_from_VFDB.tar.gz",
                 "DNA sequence of core dataset": "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz",
                 "Protein sequences of core dataset": "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz",
                 "DNA sequences of full dataset": "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz",
                 "Protein sequences of full dataset": "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz"}
    def __init__(self):
        dt = ReferenceDescriberTemplate()
        dt.name = "VFDB"
        dt.description = "A reference database for bacterial virulence factors"
        dt.documentation = "http://www.mgc.ac.cn/VFs/main.htm"
        dt.alias = "{a}_v{b}".format(a=dt.name.lower(), b=self._get_last_friday())
        self.reference_dir = "/data/reference/{a}/{b}/".format(a=dt.name, b=dt.alias).replace(" ", "_")
        dt.reference_dir = self.reference_dir
        subprocess.getoutput("rm -rf {}".format(self.reference_dir))
        os.makedirs(dt.reference_dir, exist_ok=True)
        self.raw_nfasta = "{a}{b}.fasta".format(a=dt.reference_dir, b=dt.alias)
        self.index_dir = "{}index/".format(dt.reference_dir)
        self._download()
        self._get_index_guide()
        print("Please replace the following lines in control script:")
        # refdata = "{a}{b}.json".format(a=self.index_dir, b=self.alias)
        dt.export()
        print("Please change the 'refdata' field after indexing is finished")
    @staticmethod
    def _get_last_friday():
        # VFDB updates at every Friday
        import datetime
        import calendar
        last_friday = datetime.date.today()
        oneday = datetime.timedelta(days=1)
        while last_friday.weekday() != calendar.FRIDAY:
            last_friday -= oneday
        return "{:%Y.%m.%d}".format(last_friday)
    def _download(self):
        for key in self.urls_dict:
            print("Now downloading: {}".format(key))
            url = self.urls_dict[key]
            file = url.split("/")[-1]
            cmd = "cd {A} && curl -fsSL {B} -o {C} && gunzip {C}".format(A=self.reference_dir, B=url, C=file)
            if file.endswith(".tar.gz"):
                cmd = "cd {A} && curl -fsSL {B} | tar -xz -C {A}".format(A=self.reference_dir, B=url)
            print(subprocess.getoutput(cmd))
        print(subprocess.getoutput("ln -s {A} {B}".format(A="{}VFDB_setB_nt.fas".format(self.reference_dir), B=self.raw_nfasta)))
        print("Finished downloading new raw reference: {}".format(self.raw_nfasta))
    def _get_index_guide(self):
        from meta.scripts.LaunchGuideLiner import LaunchGuideLiner
        LaunchGuideLiner.get_index_guide(index_directory=self.index_dir, raw_nfasta_file=self.raw_nfasta)


class Annotator:
    def __init__(self, annotation: str):

        self.annotation_file = annotation
        self.annotation_df = pd.read_table(self.annotation_file, sep="\t", header=0).set_index("reference_id").sort_index()
        self.annotation_df["host_strain"] = self.annotation_df["former_id"].apply(
            lambda x: re.findall("\[([^\[\]]+)\]$", x.strip())[0].strip())
        self.annotation_df["host_species"] = self.annotation_df["host_strain"].apply(
            lambda x: re.findall("([A-Z]{1}[a-z]+ [a-z]+)", x)[0].strip())
        self.annotation_df["host_genera"] = self.annotation_df["host_strain"].apply(
            lambda x: re.findall("([A-Z]{1}[a-z]+)", x)[0].strip())
        self.annotation_df["vfdb_id"] = self.annotation_df["former_id"].apply(
            lambda x: re.findall("(^VFG[0-9]+)", x)[0].strip())
        self.annotation_df["genbank_id"] = self.annotation_df["former_id"].apply(self.get_genbank_id)
        self.annotation_df["gene_name"] = self.annotation_df["former_id"].apply(self.get_gene_name)
        self.annotation_df["gene_description"] = self.annotation_df["former_id"].apply(self.get_gene_description)
        self.annotation_df = self.annotation_df.loc[:, [i for i in list(self.annotation_df) if i != "id_bp"] + ["id_bp"]]
    def update(self):
        from shutil import copy2
        copy2(self.annotation_file, "{}.bak".format(self.annotation_file))
        self.annotation_df.to_csv(self.annotation_file, sep='\t', header=True, index=True)
        print("Updated annotation: '{}'".format(self.annotation_file))
    @staticmethod
    def is_string_valid(s: str):
        s = s.strip()
        if len(s) > 0 and s != "-":
            return s
    @staticmethod
    def process_regex(pattern: str, s: str):
        out = re.findall(pattern, s)
        if len(out) > 0:
            out = Annotator.is_string_valid(out[0])
            if out:
                return out
        return np.nan
    @staticmethod
    def get_genbank_id(s: str):
        return Annotator.process_regex("^VFG[0-9]+\(([^\(\)]+)\)", s)
    @staticmethod
    def get_gene_name(s: str):
        return Annotator.process_regex("^VFG[^ ]+ \(([^\(\)]+)\)", s)
    @staticmethod
    def get_gene_description(s: str):
        return Annotator.process_regex("\)([^\[\]\(\)]+)\[", s)


if __name__ == '__main__':
    sequenceRetriever = SequenceRetriever()
    annotator = Annotator("/data/reference/VFDB/vfdb_v2018.11.09/index/vfdb_v2018.11.09_annotation.tsv")
    annotator.update()
