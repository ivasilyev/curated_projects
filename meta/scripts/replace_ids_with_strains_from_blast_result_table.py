#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from meta.utils.pandas import load_tsv
from meta.utils.io import load_string, dump_string


def parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(
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

    text_content = load_string(text_file)
    combined_blast_result_df = load_tsv(combined_blast_result_file).set_index("geninfo_id")
    print(f"Loaded replacer table with the shape {combined_blast_result_df.shape}")

    renaming_dict = combined_blast_result_df["strain"].map(lambda x: str(x).strip()).to_dict()

    text_content_replaced = re.sub("\.(gbk|gff)", "", text_content)
    counter = 0
    for renaming_key, renaming_value in renaming_dict.items():
        text_content_replaced = text_content_replaced.replace(
            *[str(i) for i in [renaming_key, renaming_value]]
        )
        counter += 1

    print(f"{counter} replacements were performed")
    dump_string(text_content_replaced, out_file)
