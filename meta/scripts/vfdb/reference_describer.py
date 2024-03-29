#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import re
import joblib as jb
import pandas as pd
from bs4 import BeautifulSoup
from time import perf_counter
from datetime import datetime
from urllib.parse import urljoin
from meta.utils.pandas import merge
from meta.utils.primitive import safe_findall
from meta.utils.language import regex_based_tokenization
from meta.utils.web import get_soup, download_file_to_dir
from meta.utils.bio_sequence import load_headers_from_fasta
from meta.utils.date_time import get_timestamp, count_elapsed_seconds
from meta.utils.file_system import decompress_file, find_file_by_tail
from meta.scripts.reference_data import AnnotatorTemplate, ReferenceData, ReferenceDescriberTemplate, SequenceRetrieverTemplate


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
    _VFDB_REGEXES = {
        "VFID": ("^([^\(\)]+)", "^([^\(\)]+)"),
        "gene_host": ("\[([^\]]+)\] *$", "(\[[^\]]+\] *$)"),
        "gene_name": ("\[([^\]]+)\] *$", "(\[[^\]]+\] *$)"),
        "gene_description": ("([^\(\)]+)$", "([^\(\)]+)$"),
        "gene_symbol": ("\(([^\(\)]+)\) *$", "^\([^\(\)]+\) *$"),
        "gene_accession_id": ("^\(([^\(\)]+)\)", "^\([^\(\)]+\) *"),
    }
    out = regex_based_tokenization(_VFDB_REGEXES, header)
    out["former_id"] = out.pop("source_string")
    out["vfdb_number"] = int(safe_findall("[0-9]+", out["VFID"]))
    return out


def mp_parse_pfasta_header(header: str):
    out = mp_parse_nfasta_header(header)
    out = {k.replace("gene", "protein"): out.get(k) for k in out}
    out["protein_header"] = out.pop("former_id")
    return out


class Annotator(AnnotatorTemplate):
    def __init__(self, refdata: ReferenceData, directory: str):
        super().__init__()
        self.refdata = refdata
        self.reference_dir = directory

        self.nucleotide_header_df = pd.DataFrame()
        self.protein_header_df = pd.DataFrame()

        self.vfs_df = pd.DataFrame()

    def load(self):
        super().load()
        start = perf_counter()
        parsed_nfasta_headers = jb.Parallel(n_jobs=-1)(
            jb.delayed(mp_parse_nfasta_header)(i) for i in self.raw_nucleotide_fasta_headers
        )
        self.nucleotide_header_df = pd.DataFrame(parsed_nfasta_headers)
        print(f"Nucleotide FASTA headers parsed into table with shape {self.nucleotide_header_df.shape} with {count_elapsed_seconds(start)}")

        start = perf_counter()
        pfasta_file = find_file_by_tail(self.reference_dir, "VFDB_setB_pro.fas")
        print(f"Use the protein FASTA file: '{pfasta_file}'")
        raw_pfasta_headers = load_headers_from_fasta(pfasta_file)
        print(f"Loaded {len(raw_pfasta_headers)} protein FASTA headers")
        parsed_pfasta_headers = jb.Parallel(n_jobs=-1)(
            jb.delayed(mp_parse_pfasta_header)(i) for i in raw_pfasta_headers
        )
        self.protein_header_df = pd.DataFrame(parsed_pfasta_headers)
        print(f"Protein FASTA headers parsed into table with shape {self.protein_header_df.shape} with {count_elapsed_seconds(start)}")

        vfs_table_file = find_file_by_tail(self.reference_dir, "VFs.xls")
        print(f"Use the VFs description file: '{vfs_table_file}'")
        self.vfs_df = pd.read_excel(vfs_table_file, engine="xlrd", header=1).fillna("")
        print(f"Loaded VFs description table with shape {self.vfs_df.shape}")
        self.vfs_df["vfdb_number"] = self.vfs_df["VFID"].str.extract("([0-9]+)").astype(int)

    def annotate(self):
        fasta_header_df = merge(self.nucleotide_header_df, self.protein_header_df, how="left",
                                on="vfdb_number")
        print(f"Merged FASTA header data into dataframe with shape {fasta_header_df.shape}")

        annotated_header_df = merge(fasta_header_df, self.vfs_df, how="left", on="vfdb_number")
        print(f"Annotated FASTA header data into dataframe with shape {annotated_header_df.shape}")

        self.annotation_df = merge(self.annotation_df, annotated_header_df, how="left",
                                   on="former_id")
        print(f"Merged final annotation dataframe with shape {self.annotation_df.shape}")


if __name__ == '__main__':
    referenceDescriber = ReferenceDescriber()
    outputDir = referenceDescriber.parse_args()
    os.makedirs(outputDir, exist_ok=True)

    sequenceRetriever = SequenceRetriever(referenceDescriber)
    sequenceRetriever.REFERENCE_ROOT_DIRECTORY = outputDir

    sequenceRetriever.get_latest_version()
    if sequenceRetriever.pick_refdata():
        print(f"Already at the latest version: '{sequenceRetriever.VERSION}'")
        startTime = perf_counter()
        annotator = Annotator(sequenceRetriever.refdata,
                              sequenceRetriever.REFERENCE_DOWNLOAD_DIRECTORY)
        annotator.run()
        annotator.validate()
        print(f"Annotation complete in {count_elapsed_seconds(startTime)}")
    else:
        print(f"Download the new version: '{sequenceRetriever.VERSION}'")
        sequenceRetriever.retrieve()
        _ = sequenceRetriever.pick_refdata()
