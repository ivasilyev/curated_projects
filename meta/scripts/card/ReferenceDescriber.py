#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} python3
"""

import os
import subprocess
import pandas as pd
from meta.templates.ReferenceDescriberTemplate import ReferenceDescriberTemplate
from meta.scripts.Utilities import Utilities


class ReferenceDescriber(ReferenceDescriberTemplate):
    NAME = "CARD"
    VERSION = "2.0.3"
    ALIAS = "card_v2.0.3"
    DESCRIPTION = "The Comprehensive Antibiotic Resistance Database"
    DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3697360/"
    WEBSITE = "https://card.mcmaster.ca/"
    REFDATA = "/data/reference/CARD/card_v2.0.3/index/card_v2.0.3_refdata.json"


class SequenceRetriever:
    _LINKS_DICT = {"Ontology Files": "https://card.mcmaster.ca/latest/ontology",
                   "Data": "https://card.mcmaster.ca/latest/data",
                   "Prevalence, Resistomes, & Variants data": "https://card.mcmaster.ca/latest/variants"}
    def __init__(self, version: str):
        self.describer = ReferenceDescriber()
        self.describer.VERSION = version
        self.describer.update_alias()
        self.reference_dir = os.path.join("/data/reference", self.describer.NAME, self.describer.ALIAS)
        self.nfasta = os.path.join(self.reference_dir, "{}.fasta".format(self.describer.ALIAS))
        self.index_dir, self.annotation_file, self._raw_nfasta_df, self._processed_nfasta_df, self.nfasta_df = (None, ) * 5
    def _download_and_unpack(self, url):
        # CARD provides data in "*.tar.bz2" archives
        url = url.strip()
        out_dir = url.strip("/").split("/")[-1].strip()
        if out_dir.count(".") > 1:
            out_dir = ".".join(out_dir.split(".")[:-1])
        out_dir = os.path.normpath(os.path.join(self.reference_dir, out_dir))
        os.makedirs(out_dir, exist_ok=True)
        cmd = "curl -fsSL {a} | tar jxf - -C {b}/".format(a=url, b=out_dir)
        print(subprocess.getoutput(cmd))
    def retrieve(self):
        if os.path.exists(self.reference_dir):
            print("Warning! The reference path exists: '{}'".format(self.reference_dir))
        os.makedirs(self.reference_dir, exist_ok=True)
        for msg in self._LINKS_DICT:
            print("Download CARD {}".format(msg))
            self._download_and_unpack(self._LINKS_DICT[msg])
        raw_nfasta = os.path.join(self.reference_dir, "data", "nucleotide_fasta_protein_homolog_model.fasta")
        if not os.path.isfile(raw_nfasta):
            raise ValueError("Not found the raw FASTA file: '{}'".format(raw_nfasta))
        print(subprocess.getoutput("ln -s {a} {b}".format(a=raw_nfasta, b=self.nfasta)))
        self.index_dir = self.describer.get_index_guide(self.nfasta)
    @staticmethod
    def _mp_parse_nfasta_header(header):
        output_dict = {"former_id": header}
        output_dict["genbank_id"] = Utilities.safe_findall("^gb\|([^|]+)", header)
        output_dict["is_antisense_strand"] = header.split("|")[2].startswith("-")
        output_dict["locus"] = Utilities.safe_findall("\|(\d+\-\d+)", header)
        output_dict["aro_id"] = Utilities.safe_findall("\|ARO:(\d+)", header)
        gene_chunk = header.split("|")[-1]
        output_dict["host"] = Utilities.safe_findall("\[(.+)\]", gene_chunk)
        output_dict["gene"] = gene_chunk.replace("[{}]".format(output_dict["host"]), "").strip()
        return Utilities.dict2pd_series(output_dict)
    def set_refdata(self, refdata_file: str):
        self.describer.REFDATA = refdata_file
    def annotate(self):
        self.annotation_file = self.describer.get_refdata_dict().get("sequence_1").annotation_file
        self._raw_nfasta_df = pd.read_table(self.annotation_file, sep='\t', header=0)
        mp_result = Utilities.multi_core_queue(self._mp_parse_nfasta_header,
                                               self._raw_nfasta_df["former_id"].values.tolist())
        self._processed_nfasta_df = Utilities.merge_pd_series_list(mp_result).sort_values("former_id")
        self.nfasta_df = Utilities.left_merge(self._raw_nfasta_df, self._processed_nfasta_df, "former_id")
        # Join 'aro_index.csv'
        aro_index_df = pd.read_table(os.path.join(self.reference_dir, "data", "aro_index.csv"), sep='\t', header=0)
        aro_index_df["aro_id"] = aro_index_df["ARO Accession"].str.extract("ARO:(\d+)")
        # 'aro_index.csv' has more entries than 'nucleotide_fasta_protein_homolog_model.fasta' provides
        self.nfasta_df = Utilities.left_merge(self.nfasta_df, aro_index_df, "aro_id")
        # Join 'aro_categories_index.csv'
        aro_categories_index_df = pd.read_table(os.path.join(self.reference_dir, "data", "aro_categories_index.csv"),
                                                sep='\t', header=0)
        self.nfasta_df = Utilities.left_merge(self.nfasta_df, aro_categories_index_df, "Protein Accession")
        # Joining 'aro_categories.csv' is useless: 'ARO Category' is filled by NaN
        # Join 'aro.csv'
        aro_df = pd.read_table(os.path.join(self.reference_dir, "ontology", "aro.csv"), sep='\t', header=0)
        aro_df.rename(columns={"Accession": "ARO Accession", "Name": "ARO Name"}, inplace=True)
        self.nfasta_df = Utilities.left_merge(self.nfasta_df, aro_df, "ARO Accession")
        self.nfasta_df = Utilities.combine_duplicate_rows(self.nfasta_df, "reference_id")
    def export_annotation(self):
        import shutil
        shutil.copy2(self.annotation_file, "{}.bak".format(self.annotation_file))
        self.nfasta_df.to_csv(self.annotation_file, sep="\t", header=True, index=False)


if __name__ == '__main__':
    # Paste updated version here
    retriever = SequenceRetriever("3.0.1")
    retriever.retrieve()
    """
    # Reference indexing (from worker node):
    
    rm -rf /data/reference/CARD/card_v3.0.1/index
    export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
    docker pull $IMG && \
    docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
    python3 /home/docker/scripts/cook_the_reference.py \
    -i /data/reference/CARD/card_v3.0.1/card_v3.0.1.fasta \
    -o /data/reference/CARD/card_v3.0.1/index
    
    # Wait until REFDATA file creates and complete the describer class template
    """
    retriever.set_refdata("/data/reference/CARD/card_v3.0.1/index/card_v3.0.1_refdata.json")
    retriever.annotate()
    retriever.export_annotation()
