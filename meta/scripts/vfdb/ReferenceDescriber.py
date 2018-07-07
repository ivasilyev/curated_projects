#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${DOCKER_IMAGE_NAME} python3
"""

"""
In the case of connection failure download data via chinese proxy from http://www.mgc.ac.cn/VFs/download.htm
VFs description file
http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz
Intra-genera VFs comparison tables
http://www.mgc.ac.cn/VFs/Down/Comparative_tables_from_VFDB.tar.gz
DNA sequence of core dataset
http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz
Protein sequences of core dataset
http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
DNA sequences of full dataset
http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
Protein sequences of full dataset
http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
"""


class ReferenceDescriber:
    name = "VFDB"
    description = "A reference database for bacterial virulence factors"
    documentation = "http://www.mgc.ac.cn/VFs/main.htm"
    # Change the following lines after reference update
    alias = "vfdb_v2018.06.29"
    refdata = "/data/reference/VFDB/vfdb_v2018.06.29/index/vfdb_v2018.06.29.json"

    def export(self):
        print("""
Database alias: {a}
REFDATA linker: {b}
              """.format(a=self.alias, b=self.refdata))


class SequenceRetriever:
    def __init__(self):
        import os
        describer = ReferenceDescriber()
        self.alias = "{a}_v{b}".format(a=describer.name.lower(), b=self._get_last_friday())
        self.reference_dir = "/data/reference/{a}/{b}/".format(a=describer.name, b=self.alias)
        os.makedirs(self.reference_dir, exist_ok=True)
        self.raw_nfasta = "{a}{b}.fasta".format(a=self.reference_dir, b=self.alias)
        self.index_dir = "{}index/".format(self.reference_dir)
        self._download()
        self._get_index_guide()
        print("Please replace the following lines in control script:")
        describer.alias = self.alias
        describer.refdata = "{a}{b}.json".format(a=self.index_dir, b=self.alias)
        describer.export()
    # VFDB updates at every Friday
    @staticmethod
    def _get_last_friday():
        import datetime
        import calendar
        last_friday = datetime.date.today()
        oneday = datetime.timedelta(days=1)
        while last_friday.weekday() != calendar.FRIDAY:
            last_friday -= oneday
        return "{:%Y.%m.%d}".format(last_friday)
    def _download(self):
        import subprocess
        url = "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz"
        cmd = "curl -fsSL {a} | zcat > {b}".format(a=url, b=self.raw_nfasta)
        print(subprocess.getoutput(cmd))
        print("Got new raw reference: {}".format(self.raw_nfasta))
    def _get_index_guide(self):
        from meta.scripts.LaunchGuideLiner import LaunchGuideLiner
        LaunchGuideLiner.get_index_guide(index_directory=self.index_dir, raw_nfasta_file=self.raw_nfasta)
