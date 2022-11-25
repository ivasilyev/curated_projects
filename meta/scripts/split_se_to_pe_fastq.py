#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Based on: https://www.biostars.org/p/19446/

import os
import re
from io import StringIO
from time import perf_counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from meta.utils.bio_sequence import dump_sequences
from meta.utils.date_time import count_elapsed_seconds


OUTPUT_FORMAT = "fastq"


def export_records(records: dict, output_mask: str):
    for strand, output_record in records.items():
        output_file = f"{output_mask}_{strand}.{OUTPUT_FORMAT}"
        dump_sequences([records[strand]], output_file, fmt=OUTPUT_FORMAT, append=True)


def split_record(record: SeqRecord):
    seq = str(record.seq)
    center = int(len(seq)/2)
    rd1_seq, rd2_seq = seq[:center], seq[center:]
    q = record.letter_annotations["phred_quality"]
    center = int(len(q)/2)
    rd1_q, rd2_q = q[:center], q[center:]
    record_id = record.id.strip("/1").strip("/2")
    rd1_id, rd2_id = f"{record_id}/1", f"{record_id}/2"
    rd1 = SeqRecord(seq=Seq(rd1_seq), id=rd1_id, description="")
    rd1.letter_annotations["phred_quality"] = rd1_q
    rd2 = SeqRecord(seq=Seq(rd2_seq), id=rd2_id, description="")
    rd2.letter_annotations["phred_quality"] = rd2_q
    return dict(R1=rd1, R2=rd2)


def read_fastq_chunk(wrapper):
    counter = 0
    buffer = ""
    while counter < 4:
        buffer += wrapper.readline()
        counter += 1
    with StringIO(buffer) as f1:
        record = list(SeqIO.parse(f1, "fastq"))[0]
        f1.close()
    if len(buffer) == 0:
        raise ValueError
    return record


def process_records(wrapper, output_mask):
    while True:
        try:
            record = read_fastq_chunk(wrapper)
            records_dict = split_record(record)
            export_records(records_dict, output_mask)
        except ValueError:
            break
    wrapper.close()


def process_sequences(file: str, output_mask: str, fmt: str = "fasta"):
    if fmt in ["fastq_gz", "fastq.gz"]:
        import mgzip
        from multiprocessing import cpu_count
        with mgzip.open(file, "rt", thread=cpu_count()) as f:
            process_records(f, output_mask)
    else:
        with open(file, mode="r", encoding="utf-8") as f:
            process_records(f, output_mask)


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

    start_0 = perf_counter()
    file_mask = os.path.join(
        output_dir,
        re.sub("(\.gz$|\.fastq\.gz$|\.fastq$|\.fq\.gz$|\.fq$)", "", os.path.basename(input_file))
    )
    process_sequences(input_file, file_mask, input_format)
    print(f"Exported SE to PE: '{input_file}' -> '{file_mask}*' in {count_elapsed_seconds(start_0)}")
