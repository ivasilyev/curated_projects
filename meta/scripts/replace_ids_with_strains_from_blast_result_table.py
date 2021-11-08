#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
from meta.scripts.Utilities import Utilities
from argparse import ArgumentParser, RawTextHelpFormatter


def parse_args():
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description="Replace GenInfo ID with strains from given table".strip(),
        epilog=""
    )
    parser.add_argument("-i", "--input", required=True,
                        help="File with text to perform replacements")
    parser.add_argument("-t", "--table", required=True,
                        help="Table with values to fetch")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file")
    _namespace = parser.parse_args()
    return _namespace.input, _namespace.table, _namespace.output


if __name__ == '__main__':
    text_file, combined_blast_result_file, out_file = parse_args()
    text_content = Utilities.load_string(text_file)
    combined_blast_result_df = Utilities.load_tsv(combined_blast_result_file).set_index(
        "geninfo_id")

    renaming_dict = combined_blast_result_df["strain"].map(
        lambda x: " ".join(Utilities.remove_empty_values(re.split("[ ]+", x)[2:]))).to_dict()

    text_content_replaced = re.sub("\.(gbk|gff)", "", text_content)
    for renaming_key, renaming_value in renaming_dict.items():
        text_content_replaced = text_content_replaced.replace(
            *(str(i) for i in [renaming_key, renaming_value]))

    text_file_data = os.path.splitext(text_file)
    Utilities.dump_string(text_content_replaced, out_file)
