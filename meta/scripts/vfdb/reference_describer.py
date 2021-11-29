#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import re
import pandas as pd
from bs4 import BeautifulSoup
from time import perf_counter
from datetime import datetime
from urllib.parse import urljoin
from meta.utils.queue import multi_core_queue
from meta.utils.primitive import safe_findall
from meta.utils.web import get_soup, download_file_to_dir
from meta.utils.bio_sequence import get_headers_from_fasta
from meta.utils.date_time import get_timestamp, count_elapsed_seconds
from meta.utils.file_system import decompress_file, find_file_by_tail
from meta.scripts.reference_data import AnnotatorTemplate, ReferenceDescriberTemplate, SequenceRetrieverTemplate


class ReferenceDescriber(ReferenceDescriberTemplate):
    NAME = "VFDB"
    DESCRIPTION = "A reference database for bacterial virulence factors"
    DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/pubmed/30395255"
    WEBSITE = "http://www.mgc.ac.cn/VFs/main.htm"


class SequenceRetriever(SequenceRetrieverTemplate):
    DOMAIN_ROOT = "http://www.mgc.ac.cn/VFs/"
    DOWNLOAD_PAGE_APPEND = "download.htm"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.download_links = ()
        self.download_page_soup = BeautifulSoup(features="lxml")

    def get_download_page_soup(self):
        download_page_soup = get_soup(urljoin(self.DOMAIN_ROOT, self.DOWNLOAD_PAGE_APPEND))
        self.download_page_soup = download_page_soup.find("table", {"id": "Table_01"})

    def get_latest_version(self):
        """
        Reference version is parsed directly from the web page.
        E.g.: 'Last update: Fri Sep 24 10:06:01 2021'
        """
        self.get_download_page_soup()
        update_text = self.download_page_soup.find("td", {"align": "right"}).find("i").text
        update_date = datetime.strptime(re.sub("^Last update: ", "", update_text), "%a %b %d %H:%M:%S %Y")
        print(f"The latest {self._reference_describer.NAME} reference version is from {get_timestamp(update_date)}")
        self.VERSION = get_timestamp(update_date, fmt="%Y.%m.%d")

    def download(self):
        for download_link in self.download_links:
            downloaded_file = download_file_to_dir(download_link, self.REFERENCE_DOWNLOAD_DIRECTORY)
            decompress_file(downloaded_file, self.REFERENCE_DOWNLOAD_DIRECTORY, remove=True)

    def retrieve(self):
        if len(self.download_page_soup) == 0:
            self.get_latest_version()
        self.download_links = [
            urljoin(self.DOMAIN_ROOT, j) for j in
            [i["href"] for i in self.download_page_soup.find_all("a", href=True)]
            if not j.endswith("htm")
        ]
        self.download()
        downloaded_nfasta = find_file_by_tail(self.REFERENCE_DOWNLOAD_DIRECTORY, "VFDB_setB_nt.fas")
        self.create_nucleotide_fasta_symlink(downloaded_nfasta, default_nfasta=True)


def mp_parse_nfasta_header(header: str):
    # Column name, regex to extract, format to remove
    _VFDB_REGEXES = {
        "vfdb_id": ("^VFG(\d+)", "VFG{}"),
        "gene_accession_id": ("\(([^\(]+)\) ", "({}) "),
        "gene_symbol": ("^\(([^\(]+)\) ", "({}) "),
        "gene_host": ("\[([^\]]+)\]$", "[{}]"),
        "gene_name": (" \[([^\]]+)\] $", " [{}] "),
        "gene_description": (".*", "{}")
    }
    out = {"former_id": header}
    # Spaces are important here
    for column_name, regexes in _VFDB_REGEXES.items():
        regex, replacement = regexes
        out[column_name] = safe_findall(regex, header, verbose=True)
        if len(out.get(column_name)) > 0:
            header = header.replace(replacement.format(out.get(column_name)), "")
    return {k: out.get(k).strip() for k in out}


def mp_parse_pfasta_header(header: str):
    out = mp_parse_nfasta_header(header)
    out = {k.replace("gene", "protein"): out.get(k) for k in out}
    out["protein_header"] = out.pop("former_id")
    return out


class Annotator(AnnotatorTemplate):
    INDEX_NAME_1 = "former_id"
    INDEX_NAME_2 = "vfdb_id"

    def __init__(self, retriever: SequenceRetriever):
        super().__init__()
        self._retriever = retriever

        self.raw_pfasta_headers = []
        self.vfs_df = pd.DataFrame()

    def load(self):
        self.refdata = self._retriever.refdata
        super().load()
        pfasta_file = find_file_by_tail(self._retriever.REFERENCE_DOWNLOAD_DIRECTORY, "VFDB_setB_pro.fas")
        print(f"Use the protein FASTA file: '{pfasta_file}'")
        self.raw_pfasta_headers = get_headers_from_fasta(pfasta_file)
        print(f"Loaded {len(self.raw_pfasta_headers)} protein FASTA headers")

        vfs_table_file = find_file_by_tail(self._retriever.REFERENCE_DOWNLOAD_DIRECTORY, "VFs.xls")
        print(f"Use the VFs description file: '{vfs_table_file}'")
        self.vfs_df = pd.read_excel(vfs_table_file, header=1).fillna("")
        print(f"Loaded VFs description file with shape '{self.vfs_df.shape}'")

    def annotate(self):
        self.load()

        raw_nfasta_headers = self.annotation_df[self.INDEX_NAME_1].values
        parsed_nfasta_headers = multi_core_queue(mp_parse_nfasta_header, raw_nfasta_headers)
        parsed_nfasta_header_df = pd.DataFrame(parsed_nfasta_headers)

        self.annotation_df = pd.concat(
            [
                i.set_index(self.INDEX_NAME_1) for i in
                [parsed_nfasta_header_df, self.annotation_df]
            ], axis=1, join="outer", sort=False
        ).rename_axis(index=self.INDEX_NAME_1).reset_index()

        parsed_pfasta_headers = multi_core_queue(mp_parse_pfasta_header, self.raw_pfasta_headers)
        parsed_pfasta_header_df = pd.DataFrame(parsed_pfasta_headers)

        zf_len = len(max(self.annotation_df[self.INDEX_NAME_2].values.tolist()))
        parsed_pfasta_header_df[self.INDEX_NAME_2] = parsed_pfasta_header_df[self.INDEX_NAME_2].str.zfill(zf_len)
        self.vfs_df[self.INDEX_NAME_2] = self.vfs_df["VFID"].str.extract("VF([0-9]+)")[0].str.zfill(zf_len)

        self.annotation_df = pd.concat(
            [
                i.set_index(self.INDEX_NAME_2).sort_index() for i in
                [self.annotation_df, parsed_pfasta_header_df, self.vfs_df]
            ], axis=1, join="outer", sort=False
        ).rename_axis(index=self.INDEX_NAME_2).sort_index().reset_index()

        self.dump()


if __name__ == '__main__':
    referenceDescriber = ReferenceDescriber()
    outputDir = referenceDescriber.parse_args()
    os.makedirs(outputDir, exist_ok=True)

    sequenceRetriever = SequenceRetriever(referenceDescriber)
    sequenceRetriever.REFERENCE_ROOT_DIRECTORY = outputDir

    sequenceRetriever.get_latest_version()
    if sequenceRetriever.pick_refdata():
        print(f"Already at the latest version: '{sequenceRetriever.VERSION}'")
        start = perf_counter()
        annotator = Annotator(sequenceRetriever)
        annotator.annotate()
        print(f"Annotation complete after {count_elapsed_seconds(start)}")
    else:
        print(f"Download new version: '{sequenceRetriever.VERSION}'")
        sequenceRetriever.retrieve()
        _ = sequenceRetriever.pick_refdata()
