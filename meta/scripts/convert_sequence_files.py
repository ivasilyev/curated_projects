#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from argparse import ArgumentParser, RawTextHelpFormatter


def parse_args():
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description="Convert sequence files from one type to another",
        epilog="Supported formats are available online: https://biopython.org/wiki/SeqIO"
    )
    parser.add_argument("-i", "--input", required=True, help="Input file")
    parser.add_argument("--input_format", required=False, default="genbank",
                        help="(Optional) Input file format")
    parser.add_argument("-o", "--output", required=True, help="Output file")
    parser.add_argument("--output_format", required=False, default="fasta",
                        help="(Optional) Output file format")
    _namespace = parser.parse_args()
    return _namespace.input, _namespace.input_format, _namespace.output, _namespace.output_format


def parse_and_write_sequences(input_file: str, input_format: str,
                              output_file: str, output_format: str):
    print(f"Converting '{input_file}' -> '{output_file}' ({input_format} -> {output_format})")
    SeqIO.write(
        SeqIO.parse(input_file, input_format),
        output_file,
        output_format
    )


if __name__ == '__main__':
    args = parse_args()
    parse_and_write_sequences(*args)
