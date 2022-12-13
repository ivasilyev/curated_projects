#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from meta.sample_data.sample_data_array import SampleDataArray

DEFAULT_READS_EXTENSION = "fastq.gz"


def parse_args():
    from argparse import ArgumentParser
    p = ArgumentParser(
        description="Generate sample data based on scanned raw reads files in the directory".strip()
    )
    p.add_argument("-i", "--input", nargs="+", help="Input directory (directories)")
    p.add_argument("-e", "--extension", default=DEFAULT_READS_EXTENSION,
                   help="Extension of reads files")
    p.add_argument("-s", "--sampledata", required=True, help="Output sample data file")
    namespace = p.parse_args()
    return (
        namespace.input,
        namespace.extension,
        namespace.sampledata,
    )


if __name__ == '__main__':
    (
        dirs,
        extension,
        sampledata
    ) = parse_args()
    arrays = dict()
    for input_dir in dirs:
        arr = SampleDataArray.import_from_dir(input_dir, extension)
        arrays[input_dir] = arr
    merged_array = SampleDataArray.merge_arrays(list(arrays.values()))
    merged_array.dump_dict(merged_array, sampledata)
