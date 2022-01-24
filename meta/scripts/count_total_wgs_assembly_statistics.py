#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from meta.scripts.Utilities import Utilities
from meta.utils.io import dump_dict
from Bio.SeqUtils import GC
from meta.utils.bio_sequence import load_sequences

reference_sequences = load_sequences("/data1/bio/projects/inicolaeva/salmonella_enterica_eclair/pga-pe-pipeline/references/blast/1002004098.gbk", fmt="genbank")
sum(len(i) for i in reference_sequences)


def parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description="Tool to count total WGS assembly statistics".strip(),
    )
    parser.add_argument("-r", "--raw_reads", required=True, nargs=2,
                        help="Raw rads")
    parser.add_argument("--raw_reads_format", default="fastq_gz",
                        help="(Optional) Raw rads format")

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

    raw_stats = Utilities.count_reads_statistics(raw_reads[0], raw_reads_format)
    assembly_stats = Utilities.count_assembly_statistics(assembly_file, assembly_format)

    reference_sequences = load_sequences(reference_file, reference_format)
    total_reference_sequence = "".join([str(i.seq) for i in reference_sequences])
    reference_stats = dict(
        total_reference_bp=len(total_reference_sequence),
        reference_gc_percentage=GC(total_reference_sequence),
    )

    coverage_stats = Utilities.count_assembly_coverages(
        raw_reads_length_sum=raw_stats["total_reads_bp"],
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
