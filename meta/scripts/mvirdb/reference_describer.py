#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from meta.scripts.reference_data import ReferenceDescriberTemplate, SequenceRetrieverTemplate
from meta.scripts.Utilities import Utilities
from argparse import ArgumentParser, RawTextHelpFormatter


class ReferenceDescriber(ReferenceDescriberTemplate):
    NAME = "MvirDB"
    DESCRIPTION = "A microbial database of protein toxins, virulence factors and antibiotic resistance genes for bio-defence applications"
    DOCUMENTATION = "https://pubmed.ncbi.nlm.nih.gov/17090593/"
    # WEBSITE = "http://mvirdb.llnl.gov/"  # Dead
    # WEBSITE = "img.jgi.doe.gov/cgi-bin/w/main.cgi"  # Dead
    WEBSITE = "https://bitbucket.org/ilya_vasilyev/mvirdb/src/master/"  # Backup


class SequenceRetriever(SequenceRetrieverTemplate):
    VERSION = "2012.04.28"
    REFERENCE_ANNOTATION = "/data/reference/MvirDB/completeMvirDBTable.txt"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


def annotate(annotation_file: str, reference_file: str):
    annotation_df = Utilities.load_tsv(annotation_file)
    _INDEX_COLUMN = "#Virulence Factor ID"

    reference_df = pd.read_csv(reference_file, engine="python", header=0,
                               sep="\t", warn_bad_lines=True, error_bad_lines=False)
    annotation_df[_INDEX_COLUMN] = annotation_df["former_id"].str.extract("^([^|]+)|").astype(int)
    annotation_df = annotation_df.merge(reference_df, how="outer", on=_INDEX_COLUMN)
    Utilities.backup_file(annotation_df)
    Utilities.dump_tsv(annotation_df, annotation_file)


def parse_args():
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter,
                            description=f"Describe and annotate {ReferenceDescriber.NAME}",
                            epilog="")
    parser.add_argument("-n", "--nfasta", metavar="<file>", required=True,
                        help="RefData file")
    parser.add_argument("-o", "--output", metavar="<dir>", required=True,
                        help="Output directory")
    namespace = parser.parse_args()
    return namespace.refdata, namespace.output


if __name__ == '__main__':
    inputNFASTAFile, outputDir = parse_args()
    os.makedirs(outputDir, exist_ok=True)

    referenceDescriber = ReferenceDescriber()
    retriever = SequenceRetriever(referenceDescriber)
    retriever.NUCLEOTIDE_FASTA = inputNFASTAFile
    retriever.REFERENCE_ROOT_DIRECTORY = outputDir
    if retriever.pick_refdata():
        annotate(retriever.REFDATA.get_sequence_dict()["annotation"], retriever.REFERENCE_ANNOTATION)
