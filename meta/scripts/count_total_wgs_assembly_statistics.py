#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
from Bio.SeqUtils import GC
from time import perf_counter
from meta.utils.io import dump_dict
from meta.utils.primitive import get_first_dict_value
from meta.utils.date_time import count_elapsed_seconds
from meta.utils.bio_sequence import load_sequences, join_sequences
from meta.scripts.Utilities import Utilities


def parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description="Tool to count total WGS assembly statistics".strip(),
        epilog="Supported formats: https://biopython.org/wiki/SeqIO"
    )
    parser.add_argument("-r", "--raw_reads", required=True, nargs="+",
                        help="Raw reads files. Any number of reads files may be supplied, only first will be used for coverage count")
    parser.add_argument("--raw_reads_format", default="fastq_gz",
                        help="(Optional) Raw reads common format")

    parser.add_argument("-a", "--assembly", required=True,
                        help="Genome assembly from these reads")
    parser.add_argument("--assembly_format", default="fasta",
                        help="(Optional) Genome assembly format")

    parser.add_argument("-e", "--reference", required=True,
                        help="Reference sequence")
    parser.add_argument("--reference_format", default="genbank",
                        help="(Optional) Reference format")

    parser.add_argument("-o", "--output", required=True,
                        help="Output file")
    _namespace = parser.parse_args()
    return (_namespace.raw_reads,
            _namespace.raw_reads_format,
            _namespace.assembly,
            _namespace.assembly_format,
            _namespace.reference,
            _namespace.reference_format,
            _namespace.output)


if __name__ == '__main__':
    (raw_reads,
     raw_reads_format,
     assembly_file,
     assembly_format,
     reference_file,
     reference_format,
     output_file) = parse_args()

    start = perf_counter()

    print(f"Counting statistics for {len(raw_reads)} raw reads files")
    raw_stats = {
        os.path.basename(i): Utilities.count_reads_statistics(i, raw_reads_format)
        for i in raw_reads
    }

    print(f"Counting statistics for WGS assembly")
    assembly_stats = Utilities.count_assembly_statistics(assembly_file, assembly_format)

    print(f"Counting statistics for reference")
    reference_sequences = load_sequences(reference_file, reference_format)
    total_reference_sequence = join_sequences(reference_sequences)
    reference_stats = dict(
        reference_file=reference_file,
        total_reference_bp=len(total_reference_sequence),
        reference_gc_percentage=GC(total_reference_sequence),
    )

    print(f"Counting coverage statistics")
    coverage_stats = Utilities.count_assembly_coverages(
        raw_reads_length_sum=get_first_dict_value(raw_stats)["total_reads_bp"],
        assembly_length=assembly_stats["total_contigs_bp"],
        reference_length=reference_stats["total_reference_bp"],
    )

    output_dict = dict(
        raw_statistics=raw_stats,
        assembly_statistics=assembly_stats,
        reference_statistics=reference_stats,
        coverage_statistics=coverage_stats,
    )

    dump_dict(output_dict, output_file)
    print(f"The calculation of WGS assembly statistics completed in {count_elapsed_seconds(start)}")
