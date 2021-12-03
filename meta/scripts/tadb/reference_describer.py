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
from meta.utils.web import get_soup, parse_table, parse_links_from_soup
from meta.utils.pandas import find_notna_rows, deduplicate_df_by_row_merging
from meta.utils.primitive import flatten_2d_array, safe_findall, remove_empty_values
from meta.scripts.reference_data import AnnotatorTemplate, ReferenceData, ReferenceDescriberTemplate, SequenceRetrieverTemplate
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

        dump_sequences(records_by_sequence_type["protein"], self.PROTEIN_FASTA, "fasta")
        print("{} nucleotide sequences were saved into the file '{}'".format(
            len(records_by_sequence_type["protein"]), self.PROTEIN_FASTA))


class Annotator(AnnotatorTemplate):
    def __init__(self, refdata: ReferenceData, nfasta_file: str, pfasta_file: str):
        super().__init__()
        self.refdata = refdata

        self._annotation_df = pd.DataFrame()

        self._nfasta_file = nfasta_file
        self._pfasta_file = pfasta_file

        self._nucleotide_headers = []
        self._protein_headers = []
        self.nucleotide_header_df = pd.DataFrame()
        self.protein_header_df = pd.DataFrame()
        self.header_based_df = pd.DataFrame()
        self.header_annotation_df = pd.DataFrame()

        self.tadb_numbers = []
        self.scrapped_feature_df = pd.DataFrame()
        self.scrapped_header_annotation_df = pd.DataFrame()

    def load(self):
        super().load()

        start = perf_counter()
        self._nucleotide_headers = load_headers_from_fasta(self._nfasta_file)
        parsed_nucleotide_headers = jb.Parallel(n_jobs=-1)(
            jb.delayed(Annotator.mp_parse_nfasta_header)
            (i) for i in self._nucleotide_headers
        )
        self.nucleotide_header_df = pd.DataFrame(parsed_nucleotide_headers).dropna(
            axis=1, how="all").dropna(axis=0, how="all").drop_duplicates()
        print(f"Parsed nucleotide header table with shape {self.nucleotide_header_df.shape}")

        self._protein_headers = load_headers_from_fasta(self._pfasta_file)
        parsed_protein_headers = jb.Parallel(n_jobs=-1)(
            jb.delayed(Annotator.mp_parse_pfasta_header)
            (i) for i in self._protein_headers
        )
        self.protein_header_df = pd.DataFrame(parsed_protein_headers).dropna(
            axis=1, how="all").dropna(axis=0, how="all").drop_duplicates()
        print(f"Parsed protein header table with shape {self.protein_header_df.shape}")
        print(f"Completed FASTA header loading in {count_elapsed_seconds(start)}")

        self.tadb_numbers = sorted(remove_empty_values(self.nucleotide_header_df["tadb_number"].unique()))
        start = perf_counter()
        tadb_scrapped_features = jb.Parallel(n_jobs=-1)(
            jb.delayed(Annotator.get_tadb_feature_table_dict)
            (i) for i in self.tadb_numbers
        )
        self.scrapped_feature_df = pd.DataFrame(tadb_scrapped_features).dropna(
            axis=1, how="all").dropna(axis=0, how="all").drop_duplicates()
        print(f"Parse web-scrapped table with shape {self.protein_header_df.shape}")
        print(f"Completed TADB feature scrapping in {count_elapsed_seconds(start)}")

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
        out["tadb_number"] = safe_findall("^[A-Z]*([0-9]+)", str(out["tadb_id"]).upper(), verbose=False)
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
        self.header_based_df = deduplicate_df_by_row_merging(
            pd.merge(
                self.nucleotide_header_df, self.protein_header_df, how="outer", on="tadb_id"
            ).drop_duplicates(),
            on="tadb_id"
        )
        print(f"Merged nucleotide & protein header data into table with shape {self.header_based_df.shape}")

        self.header_annotation_df = pd.merge(
            self.annotation_df, self.header_based_df, how="outer", on="former_id"
        )
        self.header_annotation_df = find_notna_rows(
            self.header_annotation_df, "reference_id"
        )
        print(f"Merged annotated and header data into table with shape {self.header_annotation_df.shape}")

        self.scrapped_header_annotation_df = pd.merge(
            self.header_annotation_df, self.scrapped_feature_df, how="outer", on="tadb_number"
        )
        self.scrapped_header_annotation_df = find_notna_rows(
            self.scrapped_header_annotation_df, "reference_id"
        )
        print(f"Merged annotated, header and scrapped data into table with shape {self.header_annotation_df.shape}")

        self._annotation_df = self.annotation_df.copy()
        self.annotation_df = self.scrapped_header_annotation_df


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

        sequenceRetriever.reset_nucleotide_fasta()
        sequenceRetriever.reset_protein_fasta()

        annotator = Annotator(sequenceRetriever.refdata,
                              sequenceRetriever.NUCLEOTIDE_FASTA,
                              sequenceRetriever.PROTEIN_FASTA)
        annotator.refdata = sequenceRetriever.refdata
        annotator.run()
        print(f"Annotation complete in {count_elapsed_seconds(startTime)}")
    else:
        print(f"Download the new version: '{sequenceRetriever.VERSION}'")
        sequenceRetriever.retrieve()
        _ = sequenceRetriever.pick_refdata()
