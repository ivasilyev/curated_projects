#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from time import perf_counter
from meta.utils.date_time import count_elapsed_seconds
from meta.scripts.reference_data import AnnotatorTemplate, ReferenceDescriberTemplate, SequenceRetrieverTemplate


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
    REFERENCE_ANNOTATION = "/data/reference/MvirDB/completeMvirDBTable.txt"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class Annotator(AnnotatorTemplate):
    def __init__(self, retriever: SequenceRetriever):
        super().__init__()
        self.refdata = retriever.refdata
        self.reference_annotation = retriever.REFERENCE_ANNOTATION

    def annotate(self):
        _INDEX_COLUMN = "#Virulence Factor ID"
        self.parse_annotation()
        reference_df = pd.read_csv(self.reference_annotation, engine="python", header=0,
                                   sep="\t", warn_bad_lines=True, error_bad_lines=False)
        self.annotation_df[_INDEX_COLUMN] = self.annotation_df["former_id"].str.extract("^([^|]+)|").astype(int)
        self.annotation_df = self.annotation_df.merge(reference_df, how="outer", on=_INDEX_COLUMN)
        self.dump()


if __name__ == '__main__':
    referenceDescriber = ReferenceDescriber()
    outputDir = referenceDescriber.parse_args()
    os.makedirs(outputDir, exist_ok=True)

    sequenceRetriever = SequenceRetriever(referenceDescriber)
    sequenceRetriever.REFERENCE_ROOT_DIRECTORY = outputDir
    # sequenceRetriever.retrieve()  # Not applicable
    if sequenceRetriever.pick_refdata():
        start = perf_counter()
        annotator = Annotator(sequenceRetriever)
        annotator.annotate()
        print(f"Annotation complete after {count_elapsed_seconds(count_elapsed_seconds)}")
