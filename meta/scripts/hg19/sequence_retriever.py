#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} \
bash

git pull
python3

from meta.scripts.hg19.sequence_retriever import SequenceRetriever
from meta.scripts.Utilities import Utilities
"""

import os
import pandas as pd
import gzip
from Bio import SeqIO
from .reference_describer import ReferenceDescriber
from meta.scripts.Utilities import Utilities


class SequenceRetriever:
    CHROMOSOMES = list(range(1, 23)) + ["X", "Y", "M"]

    def __init__(self, input_version: str):
        self.describer = ReferenceDescriber()
        self.describer.VERSION = input_version
        self.describer.update_alias()
        self.reference_dir = os.path.join("/data/reference/homo_sapiens/ucsc/hg", self.describer.NAME,
                                          self.describer.ALIAS)
        self.nfasta_file = os.path.join(self.reference_dir, "{}.fasta".format(self.describer.ALIAS))
        self.index_dir, self.annotation_file, self.parsed_records, self.nfasta_df = (None, ) * 4

    @staticmethod
    def _dl_wrapper(d: dict):
        url = "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr{}.fa.gz".format(d.get("chromosome"))
        return Utilities.download_file(url, d.get("chromosomes_dir"))

    @staticmethod
    def _parse_gzip_fna(gzip_fna):
        with gzip.open(gzip_fna, "rt") as f:
            seq_records = list(SeqIO.parse(f, "fasta"))
            f.close()
        return seq_records

    def retrieve(self):
        if os.path.exists(self.reference_dir):
            print("Warning! The reference path exists: '{}'".format(self.reference_dir))
        os.makedirs(self.reference_dir, exist_ok=True)
        chromosomes_dir = os.path.join(self.reference_dir, "chromosomes")
        os.makedirs(chromosomes_dir, exist_ok=True)
        # UCSC returns HTTP 530 when attempting to download in multi-thread
        compressed_chromosomes = Utilities.single_core_queue(
            self._dl_wrapper, [{"chromosome": i, "chromosomes_dir": chromosomes_dir} for i in self.CHROMOSOMES])
        # Process sequence
        self.parsed_records = Utilities.flatten_2d_array(Utilities.single_core_queue(
            self._parse_gzip_fna, compressed_chromosomes))
        self.nfasta_file = os.path.join(self.reference_dir, "hg19.fasta")
        SeqIO.write(self.parsed_records, self.nfasta_file, "fasta")
        # Process annotation
        self.index_dir = self.describer.get_index_guide(self.nfasta_file)

    def set_refdata(self, refdata_file: str):
        self.describer.set_refdata(refdata_file)

    def annotate(self):
        self.annotation_file = self.describer.get_refdata_dict().get("sequence_1").annotation_file
        annotation_dicts = [{"reference_id": i.id, "reference_bp": len(i)} for i in self.parsed_records]
        self.nfasta_df = pd.DataFrame(annotation_dicts)
        # Not required
        # Utilities.dump_tsv(annotation_df, self.annotation_file)


if __name__ == '__main__':
    # Paste updated version here
    retriever = SequenceRetriever(Utilities.get_time())
    retriever.retrieve()
    """
# Reference indexing (from worker node):
rm -rf /data/reference/homo_sapiens/ucsc/hg/hg19/hg19_v2019-07-25-14-42-32/index
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/cook_the_reference.py \
-i /data/reference/homo_sapiens/ucsc/hg/hg19/hg19_v2019-07-25-14-42-32/hg19.fasta \
-o /data/reference/homo_sapiens/ucsc/hg/hg19/hg19_v2019-07-25-14-42-32/index

# Wait until REFDATA file creates and complete the describer class template
Created reference data linker: '/data/reference/homo_sapiens/ucsc/hg/hg19/hg19_v2019-07-25-14-42-32/index/hg19_refdata.json'
    """
    retriever.set_refdata("/data/reference/homo_sapiens/ucsc/hg/hg19/hg19_v2019-07-25-14-42-32/index/hg19_refdata.json")
    """
Please update the following script lines:
class ReferenceDescriber(ReferenceDescriberTemplate):
    NAME = "hg19"
    VERSION = "2019-07-25-14-42-32"
    ALIAS = "hg19_v2019-07-25-14-42-32"
    DESCRIPTION = "The February 2009 human reference sequence (GRCh37)"
    DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/"
    WEBSITE = "https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg19&redirect=manual&source=genome.ucsc.edu"
    REFDATA = "/data/reference/homo_sapiens/ucsc/hg/hg19/hg19_v2019-07-25-14-42-32/index/hg19_refdata.json"

    """
