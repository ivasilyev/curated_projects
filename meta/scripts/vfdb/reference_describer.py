#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import re
from datetime import datetime
from meta.scripts.reference_data import AnnotatorTemplate, ReferenceDescriberTemplate, SequenceRetrieverTemplate
from meta.utils.web import get_soup, download_file_to_dir
from meta.utils.date_time import get_timestamp
from urllib.parse import urljoin
from meta.utils.file_system import decompress_file, find_file_by_tail


class ReferenceDescriber(ReferenceDescriberTemplate):
    NAME = "VFDB"
    DESCRIPTION = "A reference database for bacterial virulence factors"
    DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/pubmed/30395255"
    WEBSITE = "http://www.mgc.ac.cn/VFs/main.htm"


class SequenceRetriever(SequenceRetrieverTemplate):
    VERSION = ""
    NUCLEOTIDE_FASTA = ""
    REFERENCE_ANNOTATION = ""

    DOMAIN_ROOT = "http://www.mgc.ac.cn/VFs/"
    DOWNLOAD_PAGE_APPEND = "download.htm"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.download_links = ()

    def download(self):
        for download_link in self.download_links:
            downloaded_file = download_file_to_dir(download_link, self.REFERENCE_DOWNLOAD_DIRECTORY)
            decompress_file(downloaded_file, self.REFERENCE_DOWNLOAD_DIRECTORY, remove=True)

    def retrieve(self):
        download_page_soup = get_soup(urljoin(self.DOMAIN_ROOT, self.DOWNLOAD_PAGE_APPEND))
        download_page_table_soup = download_page_soup.find("table", {"id": "Table_01"})
        """
        Reference version is parsed directly from the web page.
        E.g.: 'Last update: Fri Sep 24 10:06:01 2021'
        """
        update_text = download_page_table_soup.find("td", {"align": "right"}).find("i").text
        update_date = datetime.strptime(re.sub("^Last update: ", "", update_text), "%a %b %d %H:%M:%S %Y")
        self.VERSION = get_timestamp(update_date)
        self.download_links = [
            urljoin(self.DOMAIN_ROOT, j) for j in
            [i["href"] for i in download_page_table_soup.find_all("a", href=True)]
            if not j.endswith("htm")
        ]
        self.download()
        downloaded_nfasta = find_file_by_tail(self.REFERENCE_DOWNLOAD_DIRECTORY, "VFDB_setB_nt.fas")
        self.create_nucleotide_fasta_symlink(downloaded_nfasta, default_nfasta=True)


class Annotator(AnnotatorTemplate):
    def __init__(self,):
        super().__init__()


if __name__ == '__main__':
    referenceDescriber = ReferenceDescriber()
    outputDir = referenceDescriber.parse_args()
    os.makedirs(outputDir, exist_ok=True)

    sequenceRetriever = SequenceRetriever(referenceDescriber)
    sequenceRetriever.REFERENCE_ROOT_DIRECTORY = outputDir
    sequenceRetriever.retrieve()
    if sequenceRetriever.pick_refdata():
        annotator = Annotator()
        annotator.annotate()
