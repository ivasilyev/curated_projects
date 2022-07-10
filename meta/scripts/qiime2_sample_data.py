#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from argparse import ArgumentParser
from meta.scripts.sample_data import create_sampledata_dict_from_dir, DEFAULT_READS_EXTENSION


def convert_sampledata(sampledata_dict: dict):
    sample_data_dicts = []
    for sampledata_line in sampledata_dict.values():
        for sampledata_reads_file, direction in zip(
            sorted(sampledata_line["reads"]), ["forward", "reverse"]
        ):
            sample_data_dicts.append({
                "sample-id": sampledata_line["name"],
                "absolute-filepath": sampledata_reads_file,
                "direction": direction
            })
    s = "categorical"
    meta_data_dicts = [{
        "#SampleID": "#q2:types",
        "BarcodeSequence": s,
        "LinkerPrimerSequence": s,
        "Description": s,
        "sample_source": s
    }, ]
    for sample_name in sorted(sampledata_dict.keys()):
        meta_data_dicts.extend([{
            "#SampleID": sample_name,
            "BarcodeSequence": "",
            "LinkerPrimerSequence": "",
            "Description": sample_name,
            "sample_source": ""
        }])
    return {
        "sample_data": pd.DataFrame(sample_data_dicts).sort_values("absolute-filepath"),
        "meta_data": pd.DataFrame(meta_data_dicts)
    }


def convert_and_dump_sampledata(sampledata_dict: dict, directory: str):
    dfs = convert_sampledata(sampledata_dict)
    os.makedirs(directory, exist_ok=True)
    for key, df in dfs.items():
        if key == "sample_data":
            sep = ","
            ext = "csv"
        else:
            sep = "\t"
            ext = "tsv"
        df.to_csv(os.path.join(directory, f"{key}.{ext}"), sep=sep, header=True, index=False)


def parse_args():
    parser = ArgumentParser(
        description="Generate sample data files for QIIME2 based on scanned raw reads files in the "
                    "input directory and dump them into the output directory".strip(),
        epilog=""
    )
    parser.add_argument("-i", "--input", nargs="+",
                        help="Input directory")
    parser.add_argument("-e", "--extension", default=DEFAULT_READS_EXTENSION,
                        help="Extension of reads files")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Output directory")
    _namespace = parser.parse_args()
    return _namespace.input, _namespace.extension, _namespace.output


if __name__ == '__main__':
    inputDir, inputExtension, outputDir = parse_args()
    sampledata_dict = create_sampledata_dict_from_dir(inputDir, inputExtension)
    convert_and_dump_sampledata(sampledata_dict, outputDir)
    print(f"Sampledata created in: '{outputDir}'")
