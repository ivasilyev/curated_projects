#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from meta.scripts.Utilities import Utilities
from argparse import ArgumentParser, RawTextHelpFormatter
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
        self.parse_annotation()

    def annotate(self):
        _INDEX_COLUMN = "#Virulence Factor ID"
        annotation_df = Utilities.load_tsv(self.annotation_file)
        reference_df = pd.read_csv(self.reference_annotation, engine="python", header=0,
                                   sep="\t", warn_bad_lines=True, error_bad_lines=False)
        annotation_df[_INDEX_COLUMN] = annotation_df["former_id"].str.extract("^([^|]+)|").astype(int)
        self.annotation_df = annotation_df.merge(reference_df, how="outer", on=_INDEX_COLUMN)
        Utilities.backup_file(self.annotation_file)
        Utilities.dump_tsv(df=self.annotation_df, table_file=self.annotation_file)


def parse_args():
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter,
                            description=f"Describe and annotate {ReferenceDescriber.NAME}",
                            epilog="")
    parser.add_argument("-o", "--output", metavar="<dir>", required=True,
                        help="Output directory")
    namespace = parser.parse_args()
    return namespace.output


if __name__ == '__main__':
    outputDir = parse_args()
    os.makedirs(outputDir, exist_ok=True)

    referenceDescriber = ReferenceDescriber()
    sequenceRetriever = SequenceRetriever(referenceDescriber)
    sequenceRetriever.REFERENCE_ROOT_DIRECTORY = outputDir
    if sequenceRetriever.pick_refdata():
        annotator = Annotator(sequenceRetriever)
        annotator.annotate()
