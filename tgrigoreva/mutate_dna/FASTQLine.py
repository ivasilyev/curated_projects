#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from tgrigoreva.mutate_dna.FASTQMutator import FASTQMutator
import re
import numpy as np


class FASTQLine:
    MUTATION_PROBABILITY_ORDER = 10
    header = ""
    sequence = ""
    header_doubler = "+"
    quality_string = ""
    quality_values = ()
    def __init__(self, header, sequence, header_doubler, quality_string):
        self.header = re.sub("^@", "", header).strip()
        self.sequence = sequence.strip().upper()
        self.header_doubler = header_doubler.strip()
        self.quality_string = quality_string.strip()
    def set_quality_values(self):
        # For the Phred33 score gauge
        # Phred quality score / ASCII table additional info:
        # http://www.somewhereville.com/2011/12/16/sanger-and-illumina-1-3-and-solexa-phred-score-q-ascii-glyph-base-error-conversion-tables/
        self.quality_values = [ord(i) - 33 for i in self.quality_string.strip()]
    def set_quality_string(self):
        self.quality_string = "".join([chr(i + 33) for i in self.quality_values])
    def randomize_header(self):
        pass
    def roll_mutation(self):
        out_sequence = ""
        out_quality_values = []
        bingo = int(self.MUTATION_PROBABILITY_ORDER / 2)
        for nucleotide, quality in zip(self.sequence, self.quality_values):
            if np.random.randint(1, self.MUTATION_PROBABILITY_ORDER + 1) != bingo or nucleotide == "N":
                out_sequence += nucleotide
                out_quality_values.append(quality)
            else:
                mutator = FASTQMutator(nucleotide, quality)
                mutated_sequence, mutated_quality = mutator.mutate()
                out_sequence += mutated_sequence
                out_quality_values.extend(mutated_quality)
        self.sequence = out_sequence
        self.quality_values = out_quality_values
    def to_str(self):
        self.set_quality_string()
        return "@{}".format("\n".join([self.header, self.sequence, self.header_doubler, self.quality_string]))
