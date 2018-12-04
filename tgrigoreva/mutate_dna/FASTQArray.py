#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from tgrigoreva.mutate_dna.FASTQLine import FASTQLine
from meta.scripts.Utilities import Utilities
import re


class FASTAArray:
    """
    This class keeps FASTQLine objects list
    """
    CHUNK_SIZE = 10 ** 6
    def __init__(self, string_to_parse: str):
        print("Initial parse of FASTQ raw string, total symbols: {}".format(len(string_to_parse)))
        string_to_parse = re.sub("[\r\n]+", "\n", string_to_parse)
        print("Split the FASTQ raw string by headers")
        self.raw_fastqs_list = string_to_parse.split("\n@")
        print("Total FASTQ reads to parse: {}".format(len(self.raw_fastqs_list)))
    def parse_fastq(self, output_file: str):
        import os
        import math
        if os.path.isfile(output_file):
            os.remove(output_file)
            print("Deleted file in order to replace it with new data: '{}'".format(output_file))
        counter = 0
        first_position = 0
        last_position = self.CHUNK_SIZE
        chunks_number = math.ceil(len(self.raw_fastqs_list) / self.CHUNK_SIZE)
        while last_position < len(self.raw_fastqs_list):
            with open(output_file, mode="a", encoding="utf-8") as f:
                f.write("{}\n".format("\n".join(Utilities.multi_core_queue(self.mp_parse_fastq_line, self.raw_fastqs_list[first_position: last_position]))))
            counter += 1
            print("Passed FASTQ parse iteration: {} (of {})".format(counter, chunks_number))
            first_position += self.CHUNK_SIZE
            last_position += self.CHUNK_SIZE
        if first_position <= len(self.raw_fastqs_list):
            with open(output_file, mode="a", encoding="utf-8") as f:
                f.write("{}\n".format("\n".join(Utilities.multi_core_queue(self.mp_parse_fastq_line, self.raw_fastqs_list[first_position: len(self.raw_fastqs_list)]))))
            print("Passed FASTQ parse last iteration")
        print("Finished parse FASTQ items: {}".format(len(self.raw_fastqs_list)))
    @staticmethod
    def mp_parse_fastq_line(single_read_fastq: str):
        fastq_list = single_read_fastq.split("\n")
        fql = FASTQLine(fastq_list[0], fastq_list[1], fastq_list[2], fastq_list[3])
        fql.set_quality_values()
        fql.roll_mutation()
        return fql.to_str()
