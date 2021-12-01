#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import re
import joblib as jb
import pandas as pd
from Bio import SeqIO
from io import StringIO
from time import perf_counter
from datetime import datetime
from urllib.parse import urljoin
from meta.utils.web import get_page
from meta.utils.pandas import dump_tsv
from meta.utils.primitive import safe_findall
from meta.utils.date_time import count_elapsed_seconds
from meta.utils.language import regex_based_tokenization
from meta.utils.web import get_soup, parse_table, parse_links_from_soup
from meta.scripts.reference_data import ReferenceDescriberTemplate, SequenceRetrieverTemplate
from meta.utils.bio_sequence import dump_sequences, load_sequences, remove_duplicate_sequences


class ReferenceDescriber(ReferenceDescriberTemplate):
    NAME = "TADB"
    DESCRIPTION = "An updated database of bacterial type II toxin-antitoxin loci"
    DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/pubmed/?term=2910666"
    WEBSITE = "http://www.mgc.ac.cn/VFs/main.htm"


def mp_parse_pfasta_header(header: str):
    out = regex_based_tokenization({
        "tadb_id": ("^TADB\|([^ ]+) *", "(^TADB\|[^ ]+ *)"),
        "protein_symbol": (" *\[([^\[\]]+)\] *$", "( *\[[^\[\]]+\] *$)"),
        "protein_host": (" *\[([^\[\]]+)\] *$", "( *\[[^\[\]]+\] *$)"),
        "protein_geninfo_id": ("^gi\|([0-9]+)\|*", "(^gi\|[0-9]+\|*)"),
        "protein_refseq_id": ("^[\| ]*ref\|([^\|]+)[\| ]*", "(^[\| ]*ref\|[^\|]+[\| ]*)"),
        "protein_description": ("(.*)", "(.*)"),
    }, header)
    out["protein_header"] = out.pop("source_string")
    return out


def mp_parse_nfasta_header(header: str):
    out = regex_based_tokenization({
        "tadb_id": ("^TADB\|([^ ]+) *", "(^TADB\|[^ ]+ *)"),
        "gene_symbol": (" *\[([^\[\]]+)\] *$", "( *\[[^\[\]]+\] *$)"),
        "gene_geninfo_id": ("^gi\|([0-9]+)[\| ]*", "(^gi\|[0-9]+[\|]*)"),
        "gene_refseq_id": ("^[\| ]*ref\|([^\|]+)[\| ]*", "(^[\| ]*ref\|[^\|]+[\| ]*)"),
        "dna_strand": ("^[\| ]*:([c]*)", "(^[\| ]*:[c]*)"),
        "start_locus": ("^([0-9]+)[ -]*", "(^[0-9]+[ -]*)"),
        "end_locus": ("^[ -]*([0-9]+)[ -]*", "(^[ -]*[0-9]+[ -]*)"),
        "gene_description": ("(.*)", "(.*)"),
    }, header)
    out["former_id"] = out.pop("source_string")
    out["is_antisense_dna_strand"] = out["dna_strand"] == "c"
    return out


def get_tadb_feature_table_dict(tadb_number: int):
    feature_page_url = urljoin(SequenceRetriever.DOMAIN_ROOT, f"feature_page.php?TAs_id={tadb_number}")
    feature_page_soup = get_soup(feature_page_url)
    feature_page_table_soup = feature_page_soup.find("table")
    if feature_page_table_soup is None or len(feature_page_table_soup) == 0:
        print(f"Cannot scrap the URL: {feature_page_url}")
        return dict()
    feature_parsed_table_dict = {f"Online Feature {k}": v[0] for k, v in parse_table(feature_page_table_soup).items()}
    feature_parsed_table_dict["tadb_number"] = feature_parsed_table_dict.pop("Online Feature TA ID")
    return feature_parsed_table_dict


class SequenceRetriever(SequenceRetrieverTemplate):
    DOMAIN_ROOT = "https://bioinfo-mml.sjtu.edu.cn/TADB2/"
    DOWNLOAD_PAGE_APPEND = "download.html"

    INDEX_COLUMN = "tadb_id"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def get_sequence_type(s: str):
        return s.split("/")[-2].strip()

    @staticmethod
    def choose_tadb_id_category(s: str):
        if s.startswith("T"):
            return "Toxin"
        if s.startswith("AT"):
            return "Antitoxin"
        if s.startswith("RE"):
            return "Regulator"

    def retrieve(self):
        """
        Reference version is parsed directly from the web page.
        E.g.:

        TADB Version: 2.0
        Last Update: June, 2017
        """

        download_page_soup = get_soup(urljoin(self.DOMAIN_ROOT, self.DOWNLOAD_PAGE_APPEND))

        update_text = download_page_soup.find("em").text
        update_datetime = datetime.strptime(re.sub("^Last Update:[ ]*", "", update_text), "%B, %Y")
        self.VERSION = update_datetime.strftime("%Y.%m")
        table_soup = download_page_soup.find("table")
        download_links = parse_links_from_soup(table_soup, self.DOMAIN_ROOT)
        sequence_types = sorted(set([self.get_sequence_type(i) for i in download_links]))

        links_by_sequence_type = {i: [j for j in download_links if self.get_sequence_type(j) == i]
                                  for i in sequence_types}
        records_by_sequence_type = {i: [] for i in sequence_types}

        for sequence_type in links_by_sequence_type.keys():
            for sequence_type_link in links_by_sequence_type[sequence_type]:
                download_buffer = get_page(sequence_type_link)
                records_by_sequence_type[sequence_type].extend(
                    list(SeqIO.parse(StringIO(download_buffer.decode("utf-8")), "fasta")))

        records_by_sequence_type = {k: remove_duplicate_sequences(v) for k, v in records_by_sequence_type.items()}

        self.reset_nucleotide_fasta()

        dump_sequences(records_by_sequence_type["nucleotide"], self.NUCLEOTIDE_FASTA)
        dump_sequences(records_by_sequence_type["protein"],
                       f"{os.path.splitext(self.NUCLEOTIDE_FASTA)[0]}.faa")


if __name__ == '__main__':
    referenceDescriber = ReferenceDescriber()
    outputDir = referenceDescriber.parse_args()
    os.makedirs(outputDir, exist_ok=True)

    sequenceRetriever = SequenceRetriever(referenceDescriber)
    sequenceRetriever.REFERENCE_ROOT_DIRECTORY = outputDir

    sequenceRetriever.get_latest_version()
    if sequenceRetriever.pick_refdata():
        pass
    else:
        print(f"Download new version: '{sequenceRetriever.VERSION}'")
        sequenceRetriever.retrieve()
        _ = sequenceRetriever.pick_refdata()
