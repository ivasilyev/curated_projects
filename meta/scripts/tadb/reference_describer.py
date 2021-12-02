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
from meta.utils.web import get_page
from meta.utils.date_time import count_elapsed_seconds
from meta.utils.language import regex_based_tokenization
from meta.utils.primitive import flatten_2d_array, safe_findall
from meta.utils.web import get_soup, parse_table, parse_links_from_soup
from meta.scripts.reference_data import AnnotatorTemplate, ReferenceDescriberTemplate, SequenceRetrieverTemplate
from meta.utils.bio_sequence import dump_sequences, load_headers_from_fasta, remove_duplicate_sequences, string_to_sequences


class ReferenceDescriber(ReferenceDescriberTemplate):
    NAME = "TADB"
    DESCRIPTION = "An updated database of bacterial type II toxin-antitoxin loci"
    DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/pubmed/?term=2910666"
    WEBSITE = "http://www.mgc.ac.cn/VFs/main.htm"


class SequenceRetriever(SequenceRetrieverTemplate):
    DOMAIN_ROOT = "https://bioinfo-mml.sjtu.edu.cn/TADB2/"
    DOWNLOAD_PAGE_APPEND = "download.html"

    INDEX_COLUMN = "tadb_id"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.links_by_sequence_type = dict()
        self.download_page_soup = BeautifulSoup(features="lxml")

    def get_download_page_soup(self):
        self.download_page_soup = get_soup(urljoin(self.DOMAIN_ROOT, self.DOWNLOAD_PAGE_APPEND))

    @staticmethod
    def get_sequence_type(s: str):
        return s.split("/")[-2].strip()

    @staticmethod
    def download_sequence(url: str):
        s = get_page(url).decode("utf-8")
        return string_to_sequences(s, "fasta")

    def get_latest_version(self):
        """
        Reference version is parsed directly from the web page.
        E.g.:

        TADB Version: 2.0
        Last Update: June, 2017
        """
        print("Looking for a new version")
        self.get_download_page_soup()
        update_text = self.download_page_soup.find("em").text
        update_datetime = datetime.strptime(re.sub("^Last Update:[ ]*", "", update_text), "%B, %Y")
        self.VERSION = update_datetime.strftime("%Y.%m")

    def retrieve(self):
        if len(self.download_page_soup) == 0:
            self.get_latest_version()
        table_soup = self.download_page_soup.find("table")
        download_links = parse_links_from_soup(table_soup, self.DOMAIN_ROOT)
        sequence_types = sorted(set([self.get_sequence_type(i) for i in download_links]))

        self.links_by_sequence_type = {
            i: sorted(set([j for j in download_links if self.get_sequence_type(j) == i]))
            for i in sequence_types
        }

        records_by_sequence_type = {i: [] for i in sequence_types}

        for sequence_type, links in self.links_by_sequence_type.items():
            print(f"Downloading {len(links)} {sequence_type} sequences")
            start = perf_counter()
            sequences = jb.Parallel(n_jobs=-1)(
                jb.delayed(SequenceRetriever.download_sequence)
                (i) for i in links
            )
            sequences = flatten_2d_array(sequences)
            records_by_sequence_type[sequence_type] = remove_duplicate_sequences(sequences)
            duplicate_number = len(sequences) - len(records_by_sequence_type[sequence_type])
            print(f"{duplicate_number} {sequence_type} sequence duplicates removed")
            print(f"Completed {sequence_type} sequence download in {count_elapsed_seconds(start)}")

        self.reset_nucleotide_fasta()

        dump_sequences(records_by_sequence_type["nucleotide"], self.NUCLEOTIDE_FASTA, "fasta")
        print("{} nucleotide sequences were saved into the file '{}'".format(
            len(records_by_sequence_type["nucleotide"]), self.NUCLEOTIDE_FASTA))

        protein_fasta = f"{os.path.splitext(self.NUCLEOTIDE_FASTA)[0]}.faa"
        dump_sequences(records_by_sequence_type["protein"], protein_fasta, "fasta")
        print("{} nucleotide sequences were saved into the file '{}'".format(
            len(records_by_sequence_type["protein"]), protein_fasta))


class Annotator(AnnotatorTemplate):
    def __init__(self, retriever: SequenceRetriever):
        super().__init__()
        self._retriever = retriever

        self.nucleotide_headers = []
        self.protein_headers = []

        self.header_based_df = pd.DataFrame()
        self.tadb_numbers = []

        self.scrapped_feature_df = pd.DataFrame()

    def load(self):
        start = perf_counter()
        nfasta_file = ""
        self.nucleotide_headers = load_headers_from_fasta(nfasta_file)
        print(f"{len(self.nucleotide_headers)} nucleotide headers were loaded")
        pfasta_file = ""
        self.protein_headers = load_headers_from_fasta(pfasta_file)
        print(f"{len(self.protein_headers)} protein headers were loaded")
        print(f"Completed FASTA header loading in {count_elapsed_seconds(start)}")

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def choose_tadb_id_category(s: str):
        if s.startswith("T"):
            return "Toxin"
        if s.startswith("AT"):
            return "Antitoxin"
        if s.startswith("RE"):
            return "Regulator"

    @staticmethod
    def _process_table_header(s: str):
        return f"Online Feature {s}"

    @staticmethod
    def get_tadb_feature_table_dict(tadb_number: int):
        feature_page_url = urljoin(SequenceRetriever.DOMAIN_ROOT,
                                   f"feature_page.php?TAs_id={tadb_number}")
        feature_page_soup = get_soup(feature_page_url)
        feature_page_table_soup = feature_page_soup.find("table")
        if feature_page_table_soup is None or len(feature_page_table_soup) == 0:
            print(f"Cannot scrap the URL: {feature_page_url}")
            return dict()
        feature_parsed_table_dict = {Annotator._process_table_header(k): v[0] for k, v in
                                     parse_table(feature_page_table_soup).items()}
        feature_parsed_table_dict["tadb_number"] = feature_parsed_table_dict.pop(
            Annotator._process_table_header("TA ID"))
        feature_parsed_table_dict[Annotator._process_table_header("Feature Replicon")] = feature_parsed_table_dict[
            Annotator._process_table_header("Feature Replicon")].replace("[Browse all TAs(s) in this replicon]").strip()
        return feature_parsed_table_dict

    def annotate(self):
        nucleotide_headers = jb.Parallel(n_jobs=-1)(
            jb.delayed(Annotator.mp_parse_nfasta_header)
            (i) for i in self.nucleotide_headers
        )
        protein_headers = jb.Parallel(n_jobs=-1)(
            jb.delayed(Annotator.mp_parse_pfasta_header)
            (i) for i in self.protein_headers
        )
        header_dataframes = [pd.DataFrame(i) for i in (nucleotide_headers, protein_headers)]
        self.header_based_df = pd.merge(*header_dataframes, how="outer", on="tadb_id",
                                        sort=False).dropna(axis=1, how="all").dropna(axis=0, how="all")
        self.header_based_df["category"] = self.header_based_df["tadb_id"].map(self.choose_tadb_id_category)
        self.header_based_df["tadb_number"] = self.header_based_df["tadb_id"].map(
            lambda x: safe_findall("^[A-Z]*([0-9]+)", str(x).upper(), verbose=False))
        self.tadb_numbers = sorted(self.header_based_df["tadb_number"].unique())

        start = perf_counter()
        tadb_scrapped_features = jb.Parallel(n_jobs=-1)(
            jb.delayed(Annotator.get_tadb_feature_table_dict)
            (i) for i in self.tadb_numbers
        )
        print(f"Completed TADB feature scrapping in {count_elapsed_seconds(start)}")

        self.scrapped_feature_df = pd.DataFrame(tadb_scrapped_features).dropna(
            axis=1, how="all").dropna(axis=0, how="all").drop_duplicates()

        self.annotation_df = pd.merge(self.header_based_df, self.scrapped_feature_df,
                                      how="outer", on="tadb_number")


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
        annotator = Annotator(sequenceRetriever)
        annotator.annotate()
        annotator.dump()
        print(f"Annotation complete in {count_elapsed_seconds(startTime)}")
    else:
        print(f"Download the new version: '{sequenceRetriever.VERSION}'")
        sequenceRetriever.retrieve()
        _ = sequenceRetriever.pick_refdata()
