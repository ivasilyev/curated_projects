#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from Bio import SeqIO
from copy import deepcopy


class ArgParser:
    def __init__(self):
        from argparse import ArgumentParser
        parser = ArgumentParser(description="Tool to merge and deduplicate chromosome and plasmid nucleotide sequences",
                                epilog="Sequences must be in FASTA format")
        parser.add_argument("-c", "--chromosome", metavar="<chromosome.fna>", required=True,
                            help="Input chromosome assembly file")
        parser.add_argument("-p", "--plasmid", metavar="<plasmid.fna>", required=True,
                            help="Input plasmid assembly file")
        parser.add_argument("-o", "--output", metavar="<output.fna>", required=True, help="Output file")
        self._namespace = parser.parse_args()
        self.valid = False
        self.chromosome = self._namespace.chromosome
        self.plasmid = self._namespace.plasmid
        self.output = self._namespace.output

    @staticmethod
    def log_and_raise(log: str):
        print(log)
        raise ValueError(log)


class Merger:
    PLASMID_MARK = " ___@PLASMID"

    def __init__(self, chromosome_file: str, plasmid_file: str):
        self.plasmid_number = 0
        self.seq_records_processed = []
        self._parse(chromosome_file)
        self._parse(plasmid_file, is_plasmid=True)
        if len(self.seq_records_processed) == 0:
            raise ValueError("Nothing to process, exit.")
        self.seq_records_processed = self.remove_duplicate_sequences(self.seq_records_processed)
        contigs_number = len(self.seq_records_processed)
        plasmid_counter = 0
        for idx, seq_record_processed in enumerate(self.seq_records_processed):
            seq_record_processed.id = "contig{}".format(str(idx + 1).zfill(len(str(contigs_number))))
            if seq_record_processed.description.endswith(self.PLASMID_MARK):
                plasmid_counter += 1
                seq_record_processed.description = "[plasmid-name=unnamed{}]".format(str(plasmid_counter).zfill(len(
                    str(self.plasmid_number))))
            else:
                seq_record_processed.description = ""

    def _parse(self, assembly_file, is_plasmid: bool = False):
        if not os.path.isfile(assembly_file):
            print("The file was not found: '{}'".format(assembly_file))
            return
        seq_records_raw = sorted(list(SeqIO.parse(assembly_file, "fasta")), key=lambda x: len(x), reverse=True)
        # NCBI does not allow to submit sequences shorter than 200 nucleotides
        seq_records_valid = [i for i in seq_records_raw if len(i) > 200]
        for idx, seq_record_raw in enumerate(seq_records_valid):
            # Example `contigs.fasta` header:
            # '>NODE_1_length_42950_cov_12.6852_component_0'
            # contig_number = int(Utilities.safe_findall("^NODE_([0-9]+)", seq_record_raw.id))
            # Processed FASTA header example:
            # >contig02 [organism=Clostridium difficile] [strain=ABDC] [plasmid-name=pABDC1] [topology=circular] [completeness=complete]
            seq_record_processed = deepcopy(seq_record_raw)
            if is_plasmid:
                self.plasmid_number += 1
                seq_record_processed.description += self.PLASMID_MARK
            self.seq_records_processed.append(seq_record_processed)

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
        SeqIO.write(self.seq_records_processed, out_file, "fasta")


if __name__ == '__main__':
    argparser = ArgParser()
    merger = Merger(argparser.chromosome, argparser.plasmid)
    merger.export(argparser.output)
