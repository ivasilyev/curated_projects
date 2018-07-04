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

http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz
http://www.mgc.ac.cn/VFs/Down/Comparative_tables_from_VFDB.tar.gz
http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz
http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
"""


class ReferenceDescriber:
    def __init__(self):
        self.database_name = ""
        self.refdata = ""
    def export(self):
        print("""
              New database name: {a}
              New REFDATA linker: {b}
              """.format(a=self.database_name, b=self.refdata))


class SequenceRetriever:
    def __init__(self):
        self.database_name = "vfdb_v{}".format(self._get_last_friday())
        self.reference_dir = "/data/reference/{}/".format(self.database_name)
        self.raw_nfasta = "{a}{b}.fasta".format(a=self.reference_dir, b=self.database_name)
        self.index_dir = "{}index/".format(self.reference_dir)
        self._download()
        self._get_index_guide()
        describer = ReferenceDescriber()
        describer.database_name = self.database_name
        describer.refdata = "{a}{b}.json".format(a=self.index_dir, b=self.database_name)
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
