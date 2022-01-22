#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import joblib as jb
import pandas as pd
from time import perf_counter
from meta.utils.pandas import merge
from meta.utils.primitive import safe_findall
from meta.utils.file_system import find_file_by_tail
from meta.utils.date_time import count_elapsed_seconds
from meta.utils.language import regex_based_tokenization
from meta.utils.bio_sequence import load_headers_from_fasta
from meta.scripts.reference_data import AnnotatorTemplate, ReferenceData, ReferenceDescriberTemplate, SequenceRetrieverTemplate


class ReferenceDescriber(ReferenceDescriberTemplate):
    NAME = "MvirDB"
    DESCRIPTION = "A microbial database of protein toxins, virulence factors and antibiotic resistance genes for bio-defence applications"
    DOCUMENTATION = "https://pubmed.ncbi.nlm.nih.gov/17090593/"
    # WEBSITE = "http://mvirdb.llnl.gov/"  # Dead
    # WEBSITE = "img.jgi.doe.gov/cgi-bin/w/main.cgi"  # Dead
    WEBSITE = "https://bitbucket.org/ilya_vasilyev/mvirdb/src/master/"  # Backup


class SequenceRetriever(SequenceRetrieverTemplate):
    VERSION = "2012.04.28"
    NUCLEOTIDE_FASTA = "/data/reference/MvirDB/mvirdb_v2012.04.28.fasta"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


def mp_parse_nfasta_header(header: str):
    # E.g. "10001|vfid|47953|vsiid|68790|ssid|SubName: Full=Leader peptidase PilD; SubName: Full=Type 4 prepilin peptidase VcpD; SubName: Full=Type IV-A prepilin peptidase PilD;"
    _REGEXES = {
        "vfid": ("^([^|]+)[ ]*\|vfid[\|]*", "^([^|]+[ ]*\|vfid[\|]*)"),
        "vsiid": ("^([^|]+)[ ]*\|vsiid[\|]*", "^([^|]+[ ]*\|vsiid[\|]*)"),
        "ssid": ("^([^|]+)[ ]*\|ssid[\|]*", "^([^|]+[ ]*\|ssid[\|]*)"),
        "feature_names": ("^[ ]*([^|]+)$", "^([ ]*[^|]+)$")
    }

    out = regex_based_tokenization(_REGEXES, header)
    out["former_id"] = out.pop("source_string")
    out.update({k: safe_findall(v, out["feature_names"]) for k, v in {
        "gene_host": "\[([^\]]+)\] *$",
        "recname_full": "[^_]*RecName:_Full=([^;]+);",
        "subname_full": "[^_]*SubName:_Full=([^;]+);",
    }.items()})
    return out


class Annotator(AnnotatorTemplate):
    def __init__(self, refdata: ReferenceData, directory: str):
        super().__init__()
        self.refdata = refdata
        self.reference_dir = directory
        self.nucleotide_header_df = pd.DataFrame()
        self.reference_df = pd.DataFrame()

    def load(self):
        super().load()
        start = perf_counter()
        nfasta_file = find_file_by_tail(self.reference_dir, "virulenceDB.nucleic.fasta")
        print(f"Use the nucleotide FASTA file: '{nfasta_file}'")
        raw_nfasta_headers = load_headers_from_fasta(nfasta_file)
        print(f"Loaded {len(raw_nfasta_headers)} nucleotide FASTA headers")
        parsed_nfasta_headers = jb.Parallel(n_jobs=-1)(
            jb.delayed(mp_parse_nfasta_header)(i) for i in raw_nfasta_headers
        )
        self.nucleotide_header_df = pd.DataFrame(parsed_nfasta_headers)
        print(f"Nucleotide FASTA headers parsed into table with shape {self.nucleotide_header_df.shape} with {count_elapsed_seconds(start)}")

        reference_file = find_file_by_tail(self.reference_dir, "completeMvirDBTable.txt")
        print(f"Use the reference description file: '{reference_file}'")
        self.reference_df = pd.read_csv(
            reference_file, engine="python", header=0, on_bad_lines="warn", sep="\t"
        ).rename(columns={"#Virulence Factor ID": "vfid"}).sort_values("vfid")
        print(f"Loaded reference description table with shape {self.reference_df.shape}")

    def annotate(self):
        annotated_header_df = merge(
            self.nucleotide_header_df, self.reference_df, on="vfid"
        )
        print(f"Annotated FASTA header data into dataframe with shape {annotated_header_df.shape}")

        self.annotation_df = merge(
            self.annotation_df, annotated_header_df, on="former_id", deduplicate=True
        )
        print(f"Merged final annotation dataframe with shape {self.annotation_df.shape}")


if __name__ == '__main__':
    referenceDescriber = ReferenceDescriber()
    outputDir = referenceDescriber.parse_args()
    os.makedirs(outputDir, exist_ok=True)

    sequenceRetriever = SequenceRetriever(referenceDescriber)
    sequenceRetriever.REFERENCE_ROOT_DIRECTORY = outputDir
    # sequenceRetriever.retrieve()  # Not applicable
    if sequenceRetriever.pick_refdata():
        annotationStart = perf_counter()
        annotator = Annotator(sequenceRetriever.refdata, sequenceRetriever.REFERENCE_ROOT_DIRECTORY)
        annotator.run()
        annotator.validate()
        print(f"Annotation complete after {count_elapsed_seconds(annotationStart)}")
