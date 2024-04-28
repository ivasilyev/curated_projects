#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from meta.sample_data.sample_data import create_sampledata_dict_from_dir, DEFAULT_READS_EXTENSION
from meta.utils.qiime import convert_and_dump_sampledata


def run(
    reads_dirs: list,
    sampledata_dir: str,
    extension: str = DEFAULT_READS_EXTENSION,
    barcode_sequence: str = "",
    linker_primer_sequence: str = ""
):
    sampledata_dict = dict()
    for reads_dir in reads_dirs:
        sampledata_dict.update(create_sampledata_dict_from_dir(reads_dir, extension))
    sampledata_files = convert_and_dump_sampledata(
        directory=sampledata_dir,
        sample_data_dict=sampledata_dict,
        barcode_sequence=barcode_sequence,
        linker_primer_sequence=linker_primer_sequence,
    )
    return dict(
        sampledata_dict=sampledata_dict,
        sampledata_files=sampledata_files
    )


def parse_args():
    parser = ArgumentParser(
        description="Generate sample data files for QIIME2 based on scanned raw reads files in the "
                    "input directory and dump them into the output directory".strip(),
        epilog=""
    )
    parser.add_argument("-i", "--input", nargs="+", help="Input directories")
    parser.add_argument("-e", "--extension", default=DEFAULT_READS_EXTENSION,
                        help="Extension of reads files")
    parser.add_argument("-b", "--barcode", default="", help="Barcode sequence")
    parser.add_argument("-l", "--linker", default="", help="Linker primer sequence")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    _namespace = parser.parse_args()
    return (
        _namespace.input,
        _namespace.extension,
        _namespace.barcode,
        _namespace.linker,
        _namespace.output,
    )


if __name__ == '__main__':
    (
        inputDirs,
        inputExtension,
        inputBarcode,
        inputLinker,
        outputDir
    ) = parse_args()
    run(
        reads_dirs=inputDirs,
        extension=inputExtension,
        sampledata_dir=outputDir,
        barcode_sequence=inputBarcode,
        linker_primer_sequence=inputLinker
    )
