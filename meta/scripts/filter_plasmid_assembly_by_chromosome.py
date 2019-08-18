#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO


class ArgParser:
    def __init__(self):
        from argparse import ArgumentParser
        parser = ArgumentParser(description="Tool to remove from the plasmid nucleotide assembly the sequences which "
                                            "are also presented in chromosome assembly",
                                epilog="Sequences must be in FASTA format")
        parser.add_argument("-c", "--chromosome", metavar="<chromosome.fna>", required=True,
                            help="Input chromosome assembly file")
        parser.add_argument("-p", "--plasmid", metavar="<plasmid.fna>", required=True,
                            help="Input plasmid assembly file")
        parser.add_argument("-r", "--rename", default=False, action="store_true",
                            help="If checked, the output plasmid contigs will be renamed as 'contig<number>'")
        parser.add_argument("-o", "--output", metavar="<output.fna>", required=True,
                            help="Output plasmid assembly file")
        self._namespace = parser.parse_args()
        self.chromosome = self._namespace.chromosome
        self.plasmid = self._namespace.plasmid
        self.rename = self._namespace.rename
        self.output = self._namespace.output


class Filter:
    def __init__(self, chromosome_file: str, plasmid_file: str, rename: bool = True):
        chromosome_sequences = [i.seq for i in self._parse(chromosome_file)]
        plasmid_records = self._parse(plasmid_file)
        self.processed_records = []
        for plasmid_record in plasmid_records:
            if plasmid_record.seq not in chromosome_sequences:
                self.processed_records.append(plasmid_record)
        self.processed_records = self.remove_duplicate_sequences(self.processed_records)
        if rename:
            for idx, seq_record_processed in enumerate(self.processed_records):
                seq_record_processed.id = "contig{}".format(str(idx + 1).zfill(len(str(len(self.processed_records)))))
                seq_record_processed.description = ""

    @staticmethod
    def _parse(assembly_file):
        return sorted(list(SeqIO.parse(assembly_file, "fasta")), key=lambda x: len(x), reverse=True)

    @staticmethod
    def remove_duplicate_sequences(records: list):
        out = []
        sequences = []
        for record in records:
            if record.seq not in sequences:
                sequences.append(record.seq)
                out.append(record)
        return out

    def export(self, out_file):
        import os
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        SeqIO.write(self.processed_records, out_file, "fasta")


if __name__ == '__main__':
    argparser = ArgParser()
    filter_ = Filter(argparser.chromosome, argparser.plasmid, argparser.rename)
    filter_.export(argparser.output)
