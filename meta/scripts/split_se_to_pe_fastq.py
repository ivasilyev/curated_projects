#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Based on: https://www.biostars.org/p/19446/

import os
import re
import joblib as jb
from time import perf_counter
from Bio.SeqRecord import SeqRecord
from meta.utils.date_time import count_elapsed_seconds
from meta.utils.bio_sequence import load_sequences, dump_sequences


def process_record(record):
    seq = str(record.seq)
    rd1_seq = seq[:len(seq)/2]
    rd2_seq = seq[len(seq)/2:]
    Q = record.letter_annotations["phred_quality"]
    rd1_q = Q[:len(Q)/2]
    rd2_q = Q[len(Q)/2:]
    rd1_id = record.id.strip("/1").strip("/2") + "/1"
    rd2_id = record.id.strip("/1").strip("/2") + "/2"
    rd1 = SeqRecord(rd1_seq, id=rd1_id, description="")
    rd1.letter_annotations["phred_quality"] = rd1_q
    rd2 = SeqRecord(rd2_seq, id=rd2_id, description="")
    rd2.letter_annotations["phred_quality"] = rd2_q

    return dict(R1=rd1.format("fastq"), R2=rd2.format("fastq"))


def _parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="")
    parser.add_argument("-i", "--input_file", metavar="<file>", required=True,
                        help="Input file")
    parser.add_argument("-f", "--format", metavar="<str>", default="fastq_gz",
                        help="Sequence format")
    parser.add_argument("-o", "--output_dir", metavar="<dir>", required=True,
                        help="Output directory")
    _namespace = parser.parse_args()
    return (
        _namespace.input_file,
        _namespace.format,
        _namespace.output_dir,
    )


if __name__ == '__main__':
    (
        input_file,
        input_format,
        output_dir,
    ) = _parse_args()

    output_format = "fastq"
    start_0 = perf_counter()
    records = load_sequences(input_file, fmt=input_format, is_filter=False, is_sort=False)
    print(f"Parsed {len(records)} sequences from '{input_file}' in {count_elapsed_seconds(start_0)}")

    start_1 = perf_counter()
    file_mask = re.sub("(\.gz$|\.fastq\.gz$|\.fastq$|\.fq\.gz$|\.fq$)", "", os.path.basename(input_file))
    output_files = {f"R{i}": os.path.join(output_dir, f"{file_mask}_R{i}.{output_format}") for i in [1, 2]}
    processed_record_dicts = jb.Parallel(n_jobs=-1)(jb.delayed(process_record)(i) for i in records)
    print(f"Split SE to PE in {count_elapsed_seconds(start_1)}")

    start_2 = perf_counter()
    for strand, output_file in output_files.values():
        processed_records = [i[strand] for i in processed_record_dicts]
        dump_sequences(processed_records, output_file, fmt=output_format)
    print(f"Exported SE to PE: '{input_file}' -> '{output_files}' in {count_elapsed_seconds(start_2)}")

    print(f"Done in {count_elapsed_seconds(start_0)}")
