#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NCBI contamination report example:

Contamination_Kleb102_genome.txt

SUBID     	BioProject	BioSample	Organism
--------------------------------------------------------
SUB5915929	PRJNA556398	SAMN12349658	Klebsiella pneumoniae

[] We ran your sequences through our Contamination Screen. The screen found
contigs that need to be trimmed and/or excluded. Please adjust the
sequences appropriately and then resubmit your sequences. After you remove the
contamination, trim any Ns at the ends of the sequence and remove any sequences
that are shorter than 200 nt and not part of a multi-component scaffold.

Note that hits in eukaryotic genomes to mitochondrial sequences can be ignored
when specific criteria are met. Those criteria are explained below.

Note that mismatches between the name of the adaptor/primer identified in the screen
and the sequencing technology used to generate the sequencing data should not be used
to discount the validity of the screen results as the adaptors/primers of many
different sequencing platforms share sequence similarity.



Screened 106 sequences, 5,297,684 bp.
2 sequences to exclude, 3 sequences with locations to mask/trim

Exclude:
Sequence name, length, apparent source
NODE_51_length_527_cov_1.2	527	Pan troglodytes
NODE_85_length_443_cov_1.46519	443	Primates


Trim:
Sequence name, length, span(s), apparent source
NODE_25_length_2931_cov_1.15514	2931	1..264	Homo sapiens
NODE_40_length_793_cov_0.833333	793	1..139	Pan troglodytes
NODE_61_length_493_cov_0.770492	493	1..132	Pan troglodytes


Duplicated:
Sequence names, length
lcl|contig099 lcl|contig153 (452 bp)
lcl|contig053 lcl|contig136 (2978 bp)
contig18 RC(contig45) (565 bp)

"""
# Duplicated entry was taken from another file

from time import perf_counter
from meta.utils.date_time import count_elapsed_seconds
from meta.scripts.Utilities import Utilities
from Bio import SeqIO


class ArgParser:
    def __init__(self):
        from argparse import ArgumentParser
        parser = ArgumentParser(description="Tool to remove contamination from given FASTA nucleotide sequence file "
                                            "based on NCBI contamination report",
                                epilog="Note if the FASTA headers contain spaces, only the first piece of split header "
                                       "will be processed")
        parser.add_argument("-i", "--input", metavar="<input.fna>", required=True, help="Input FASTA file")
        parser.add_argument("-c", "--contamination", metavar="<Contamination_input.txt>", required=True,
                            help="NCBI contamination report file")
        parser.add_argument("-o", "--output", metavar="<output.fna>", required=True, help="Output FASTA file")
        self._namespace = parser.parse_args()
        self.input = self._namespace.input
        self.contamination = self._namespace.contamination
        self.output = self._namespace.output


class ContaminationRemover:
    _NCBI_MIN_SEQ_LENGTH = 200

    def __init__(self, fna_file, contamination_report):
        self.sequence_file = fna_file
        self.contamination_file = contamination_report
        self.contamination_lines = Utilities.load_2d_array(self.contamination_file)

        exclude_index = self.find_index(self.contamination_lines, ["Exclude:"])
        trim_index = self.find_index(self.contamination_lines, ["Trim:"])
        duplicated_index = self.find_index(self.contamination_lines, ["Duplicated:"])
        # Issue lines order: exclude, trim, duplicated
        headers_to_remove = []
        if exclude_index:
            exclude_lines = []
            if trim_index:
                exclude_lines = self.contamination_lines[exclude_index + 2: trim_index]
            elif duplicated_index:
                exclude_lines = self.contamination_lines[exclude_index + 2: duplicated_index]
            else:
                exclude_lines = self.contamination_lines[exclude_index + 2:]
            exclude_lines = Utilities.remove_empty_values([i[0] for i in exclude_lines])
            headers_to_remove.extend(exclude_lines)

        if trim_index:
            trim_lines_processed = dict()
            if duplicated_index:
                trim_lines_raw = self.contamination_lines[trim_index + 2: duplicated_index]
            else:
                trim_lines_raw = self.contamination_lines[trim_index + 2:]
            for trim_line_raw in Utilities.remove_empty_values(trim_lines_raw):
                for trim_span in Utilities.remove_empty_values(trim_line_raw[2].split(",")):
                    trim_indices = [int(i.strip()) for i in trim_span.split("..")]
                    # It seems that reported sequence positions are not zero-based
                    trim_indices[0] -= 1
                    trim_lines_processed[trim_line_raw[0]] = trim_indices
            headers_to_remove.extend(list(trim_lines_processed.keys()))

        if duplicated_index:
            processed_duplicated_lines = list()
            duplicated_lines = [i for i in self.contamination_lines[duplicated_index + 2:]]
            # Removing only the first occurrence in Duplicates
            for duplicated_line in duplicated_lines:
                duplicated_str = duplicated_line[0]
                if duplicated_str.startswith("# "):
                    continue
                processed_duplicated_lines.append(duplicated_str.strip().split(" ")[0])
            headers_to_remove.extend(processed_duplicated_lines)

        self.seq_records = list(SeqIO.parse(self.sequence_file, "fasta"))
        print(f"Imported {len(self.seq_records)} raw records from '{self.sequence_file}'")
        headers_to_remove = sorted(set(Utilities.remove_empty_values(headers_to_remove)))
        print("{} headers were marked to remove: '{}'".format(
            len(headers_to_remove), "', '".join(headers_to_remove)
        ))
        out_records = []
        for record_raw in Utilities.remove_duplicate_sequences(self.seq_records):
            record_id = record_raw.id.split(" ")[0].strip()
            if record_id not in headers_to_remove:
                out_records.append(record_raw)
            # Some positions from the "Trim" entry after NCBI processing were moved into the "Exclude" entry
            # The point is removing them instead of trimming
        self.valid_records = [i for i in out_records if len(i) >= self._NCBI_MIN_SEQ_LENGTH]

    @staticmethod
    def find_index(items: list, item):
        if item in items:
            return items.index(item)
        return None

    def export(self, output_file):
        SeqIO.write(self.valid_records, output_file, "fasta")
        print(f"Exported {len(self.valid_records)} valid records into '{output_file}'")


if __name__ == '__main__':
    mainParser = ArgParser()
    start = perf_counter()
    remover = ContaminationRemover(mainParser.input, mainParser.contamination)
    remover.export(mainParser.output)
    print(f"Decontaminated in {count_elapsed_seconds(start)}")
