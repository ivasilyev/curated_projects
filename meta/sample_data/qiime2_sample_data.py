#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from argparse import ArgumentParser
from meta.sample_data.sample_data import create_sampledata_dict_from_dir, DEFAULT_READS_EXTENSION


_Q2_CAT_TYPE = "categorical"


def convert_sampledata(
    sample_data_dict: dict,
    barcode_sequence: str = "",
    linker_primer_sequence: str = "",
):
    sample_data_dicts = []
    for sampledata_line in sample_data_dict.values():
        for sampledata_reads_file, direction in zip(
            sorted(sampledata_line["reads"]), ["forward", "reverse"]
        ):
            sample_data_dicts.append({
                "sample-id": sampledata_line["name"],
                "absolute-filepath": sampledata_reads_file,
                "direction": direction
            })
    meta_data_dicts = [{
        "#SampleID": "#q2:types",
        "BarcodeSequence": _Q2_CAT_TYPE,
        "LinkerPrimerSequence": _Q2_CAT_TYPE,
        "Description": _Q2_CAT_TYPE,
        "SampleSource": _Q2_CAT_TYPE
    }, ]
    for sample_name in sorted(sample_data_dict.keys()):
        meta_data_dicts.extend([{
            "#SampleID": sample_name,
            "BarcodeSequence": barcode_sequence,
            "LinkerPrimerSequence": linker_primer_sequence,
            "Description": sample_name,
            "SampleSource": ""
        }])
    return {
        "sample": pd.DataFrame(sample_data_dicts).sort_values("absolute-filepath"),
        "meta": pd.DataFrame(meta_data_dicts)
    }


def convert_and_dump_sampledata(directory: str, *args, **kwargs):
    dfs = convert_sampledata(*args, **kwargs)
    os.makedirs(directory, exist_ok=True)
    out = dict()
    for key, df in dfs.items():
        if key == "sample":
            sep = ","
            ext = "csv"
        else:
            sep = "\t"
            ext = "tsv"
        file = os.path.join(directory, f"qiime2_{key}_data.{ext}")
        df.to_csv(
            file,
            sep=sep,
            header=True,
            index=False
        )
        out[key] = file
    print(f"Sampledata created in: '{directory}'")
    return out


def run(
    reads_dirs: list,
    extension: str,
    sampledata_dir: str,
    barcode_sequence: str,
    linker_primer_sequence: str
):
    sampledata_dict = dict()
    for reads_dir in reads_dirs:
        sampledata_dict.update(create_sampledata_dict_from_dir(reads_dir, extension))
    convert_and_dump_sampledata(
        directory=sampledata_dir,
        sample_data_dict=sampledata_dict,
        barcode_sequence=barcode_sequence,
        linker_primer_sequence=linker_primer_sequence,
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
        reads_dir=inputDirs,
        extension=inputExtension,
        sampledata_dir=outputDir,
        barcode_sequence=inputBarcode,
        linker_primer_sequence=inputLinker
    )
