#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import pandas as pd

"""
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3

"""

"""
Links:

1. Sequence data for gene catalogs
Integrated non-redundant gene catalog (IGC, nucleotide sequences, fasta)
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.fa.gz

Integrated non-redundant gene catalog (IGC, amino acid sequences, fasta)
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.pep.gz

2. Gene profile, genus profile, KO profile
Gene abundance profile table for 1267 samples (9.9 million genes X 1267 samples)
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/2.ProfileTables/1267sample.gene.relativeAbun.table.gz

Genus profile table for 1267 samples
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/2.ProfileTables/IGC.genus.normalization.ProfileTable.gz

KO profile table for 1267 samples
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/2.ProfileTables/IGC.KO.normalization.ProfileTable.gz

3. Gene annotation summary
IGC annotation and occurrence frequency summary table
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/3.IGC.AnnotationInfo/IGC.annotation_OF.summary.gz

MetaHIT 2010 KEGG annotation
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/6.PreviousGeneCatalog/PGC.gene.KO.list.gz

"""


class ReferenceDescriber:
    name = "CARD"
    description = "The Comprehensive Antibiotic Resistance Database"
    documentation = "https://card.mcmaster.ca/"
    # Change the following lines after reference update
    alias = ""
    refdata = ""
    def export(self):
        print("""
Database alias: {a}
REFDATA linker: {b}
              """.format(a=self.alias, b=self.refdata))
    def parse_refdata(self):
        from meta.scripts.RefDataParser import RefDataParser
        return RefDataParser(self.refdata).get_parsed_list()


class SequenceRetriever:
    def __init__(self):
        describer = ReferenceDescriber()
        self.alias = "card_v2.0.3"
        self.reference_dir = "/data/reference/{a}/{b}/".format(a=describer.name, b=self.alias)
        os.makedirs(self.reference_dir, exist_ok=True)
        self._dl_dict = {"Ontology Files": "https://card.mcmaster.ca/latest/ontology",
                         "Data": "https://card.mcmaster.ca/latest/data",
                         "Prevalence, Resistomes, & Variants data": "https://card.mcmaster.ca/latest/variants"}
        for msg in self._dl_dict:
            print("Download CARD {}".format(msg))
            self._download_and_unpack(self._dl_dict[msg])
        self.raw_nfasta = "{}data/nucleotide_fasta_protein_homolog_model.fasta".format(self.reference_dir)
        self.processed_nfasta = "{a}{b}.fasta".format(a=self.reference_dir, b=self.alias)
        print(subprocess.getoutput("ln -s {a} {b}").format(a=self.raw_nfasta, b=self.processed_nfasta))
        self.index_dir = "{}index/".format(self.reference_dir)
        self._get_index_guide()
        print("Please replace the following lines in control script:")
        describer.alias = self.alias
        describer.refdata = "{a}{b}.json".format(a=self.index_dir, b=self.alias)
        describer.export()
    def _download_and_unpack(self, url):
        url = url.strip()
        out_dir = url.split("/")[-1]
        if out_dir.count(".") > 1:
            out_dir = ".".join(out_dir.split(".")[:-1])
        out_dir = self.reference_dir + out_dir + "/"
        os.makedirs(out_dir, exist_ok=True)
        cmd = "curl -fsSL {a} | tar jxf - -C {b}".format(a=url, b=out_dir)
        print(subprocess.getoutput(cmd))
    def _get_index_guide(self):
        from meta.scripts.LaunchGuideLiner import LaunchGuideLiner
        LaunchGuideLiner.get_index_guide(index_directory=self.index_dir, raw_nfasta_file=self.processed_nfasta)


if __name__ == '__main__':
    retriever = SequenceRetriever()


# self_reference_dir = "/data/reference/CARD/test/"
# os.makedirs(self_reference_dir, exist_ok=True)
#
# url = "https://card.mcmaster.ca/latest/data"
# url = url.strip()
# out_dir = url.strip().split("/")[-1]
# if out_dir.count(".") > 1:
#     out_dir = ".".join(out_dir.split(".")[:-1])
#
# out_dir = self_reference_dir + out_dir + "/"
# os.makedirs(out_dir, exist_ok=True)
# cmd = "curl -fsSL {a} | tar jxf - -C {b}".format(a=url, b=out_dir)
# subprocess.getoutput(cmd)
# self_raw_nfasta = "{}data/nucleotide_fasta_protein_homolog_model.fasta".format(self_reference_dir)
# self_index_dir = "{}index/".format(self_reference_dir)


"""
# Reference indexing (from worker node):

rm -rf /data/reference/IGC/igc_v2014.03/index
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && docker pull $IMG && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG python3 /home/docker/scripts/cook_the_reference.py -i /data/reference/IGC/igc_v2014.03.fasta -n -o /data/reference/IGC/igc_v2014.03/index

Wait until REFDATA file creates

Please replace the following lines in control script:

Database alias: igc_v2014.03
REFDATA linker: /data/reference/IGC/igc_v2014.03/index/igc_v2014.03_refdata.json

"""
