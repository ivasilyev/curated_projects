#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


class FASTQMutator:
    """
    All instance methods always return the tuple (<nucleotide>, [<quality>])
    """
    MAX_INSERTION_LENGTH = 5
    nucleotide = ""
    quality = 0
    def __init__(self, nucleotide: str, quality: int):
        self.nucleotide = nucleotide.strip().upper()
        self.quality = quality
    def mutate(self):
        scenario_id = np.random.randint(0, 5)
        if scenario_id < 3:
            return self._make_snp()
        elif scenario_id == 3:
            return self._make_insertion()
        else:
            return self._make_deletion()
    def _make_snp(self):
        out = self.get_random_nucleotide(self.nucleotide, self.quality)
        return out[0], [out[1]]
    def _make_insertion(self):
        insertion_length = np.random.randint(0, self.MAX_INSERTION_LENGTH + 1)
        insertion_list = [self.get_random_nucleotide(quality=self.quality) for i in np.arange(0, insertion_length)]
        nucleotide_out = "".join([i[0] for i in insertion_list])
        quality_out = [i[1] for i in insertion_list]
        return nucleotide_out, quality_out
    def _make_deletion(self):
        return "", []
    @staticmethod
    def get_random_nucleotide(nucleotide: str = "", quality: int = None):
        if not quality:
            # 93 seems to be a maximum of the Phred gauge. High integer is exclusive for np.random.randint()
            quality = np.random.randint(0, 94)
        return (np.random.choice([i for i in ("A", "T", "G", "C") if i != nucleotide], replace=False),
                abs(np.random.randint(quality - 10, quality + 11)))
