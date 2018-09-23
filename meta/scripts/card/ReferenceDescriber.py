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


class ReferenceDescriber:
    name = "CARD"
    description = "The Comprehensive Antibiotic Resistance Database"
    documentation = "https://card.mcmaster.ca/"
    # Change the following lines after reference update
    alias = "card_v2.0.3"
    refdata = "/data/reference/CARD/card_v2.0.3/index/card_v2.0.3_refdata.json"
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
        print(subprocess.getoutput("ln -s {a} {b}".format(a=self.raw_nfasta, b=self.processed_nfasta)))
        self.index_dir = "{}index/".format(self.reference_dir)
        self._get_index_guide()
        print("Please update the following lines in describer:")
        describer.alias = self.alias
        describer.refdata = "{a}{b}_refdata.json".format(a=self.index_dir, b=self.alias)
        self.annotation = "{a}{b}_annotation.tsv".format(a=self.index_dir, b=self.alias)
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
    def annotate(self):
        annotation_df = pd.read_table(self.annotation, sep='\t', header=0)
        annotation_df["aro"] = annotation_df["former_id"].str.extract("\|ARO:([0-9]+)\|")
        annotation_df["gb"] = annotation_df["former_id"].str.extract("^gb\|([^|]+)\|")
        # reference_annotation_aro_categories_df = pd.read_table("/data/reference/CARD/card_v2.0.3/data/aro_categories.csv", sep='\t', header=0)
        # reference_annotation_aro_categories_df["aro"] = reference_annotation_aro_categories_df["ARO Accession"].str.extract("ARO:([0-9]+)")
        # annotation_df = annotation_df.merge(reference_annotation_aro_categories_df, how="left", on="aro")
        # Empty merge result
        reference_annotation_index_df = pd.read_table("{}data/aro_index.csv".format(self.reference_dir), sep='\t', header=0)
        reference_annotation_index_df["aro"] = reference_annotation_index_df["ARO Accession"].str.extract("ARO:([0-9]+)")
        annotation_df = annotation_df.merge(reference_annotation_index_df, how="left", on="aro")
        #
        reference_annotation_categories_index_df = pd.read_table("{}data/aro_categories_index.csv".format(self.reference_dir), sep='\t', header=0)
        annotation_df = annotation_df.merge(reference_annotation_categories_index_df.loc[:, ["Protein Accession"] + [i for i in list(reference_annotation_categories_index_df) if i not in list(annotation_df)]], how="left", on="Protein Accession")
        annotation_df.to_csv(self.annotation, sep='\t', header=True, index=False)


if __name__ == '__main__':
    retriever = SequenceRetriever()
    retriever.annotate()
