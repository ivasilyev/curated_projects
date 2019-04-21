#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${DOCKER_IMAGE_NAME} python3
"""

import os
import re
import subprocess
import pandas as pd
from meta.templates.ReferenceDescriberTemplate import ReferenceDescriberTemplate
from meta.scripts.Utilities import Utilities


class ReferenceDescriber(ReferenceDescriberTemplate):
    NAME = "TADB"
    VERSION = "2.0"
    ALIAS = "tadb_v2.0"
    DESCRIPTION = "An updated database of bacterial type II toxin-antitoxin loci"
    DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/pubmed/?term=29106666"
    WEBSITE = "http://202.120.12.135/TADB2/index.php"
    REFDATA = "/data/reference/TADB/tadb_v2.0/index/tadb_v2.0_refdata.json"


class SequenceRetriever:
    def __init__(self, version: str):
        self.describer = ReferenceDescriber()
        self.describer.VERSION = version
        self.describer.update_alias()
        self.reference_dir = os.path.join("/data/reference", self.describer.NAME, self.describer.ALIAS)
        links = [i for i in self.scrap_links_from_web_page("http://202.120.12.135/TADB2/download.html") if
                 i.endswith(".fas")]
        self._fasta_types = ["nucleotide", "protein"]
        self._links_dict = {k: [i for i in links if i.split("/")[-2] == k] for k in self._fasta_types}
        assert len(self._links_dict["nucleotide"]) + len(self._links_dict["protein"]) == len(links)
        self.nfasta = os.path.join(self.reference_dir, "{}.fasta".format(self.describer.ALIAS))
        self.pfasta = os.path.join(self.reference_dir, "{}_protein.fasta".format(self.describer.ALIAS))
        self.index_dir = ""
    @staticmethod
    def scrap_links_from_web_page(url: str) -> list:
        import bs4
        import lxml
        soup = bs4.BeautifulSoup(Utilities.get_page(url), "lxml")
        root = "/".join(url.split("/")[:-1])
        # Scrap the web page roughly and get all required links
        raw_links = [a["href"].strip() for a in soup.find_all("a", href=True)]
        return [i if "://" in i else "{}/{}".format(root, i) for i in raw_links]
    @staticmethod
    def download_file(d: dict):
        """
        :param d: {"url": str, "out_dir": str}
        :return: None
        """
        from time import sleep
        _RETRIES_LEFT = 5
        _SLEEP_SECONDS = 3
        url = d.get("url").strip()
        out_dir = os.path.normpath(d.get("out_dir").strip())
        assert len(url) > 0 and len(out_dir) > 0
        os.makedirs(out_dir, exist_ok=True)
        while _RETRIES_LEFT > 0:
            out_file = os.path.join(out_dir, url.split("/")[-1])
            print(subprocess.getoutput("curl -fsSL {} -o {}".format(url, out_file)))
            sleep(_SLEEP_SECONDS)
            if os.path.isfile(out_file):
                print("Download finished: '{}'".format(out_file))
                break
            _RETRIES_LEFT -= 1
            print("Warning! Failed download: '{}'. Retries left: {}".format(url, _RETRIES_LEFT))
    @staticmethod
    def ls(dir_name: str):
        out = []
        for root, dirs, files in os.walk(dir_name):
            for file in files:
                out.append(os.path.join(root, file))
        return out
    @staticmethod
    def cat(*args):
        """
        :param args: concatenating files and tailing target file
        :return: None
        """
        if len(args) < 2:
            raise ValueError("Not enough files: '{}'".format(args))
        sources = sorted(args[:-1])
        with open(args[-1], mode="w", encoding="utf-8") as target:
            for source_file in sources:
                with open(source_file, mode="r", encoding="utf-8") as source:
                    for line in source:
                        target.write(line)
                    target.write("\n")
                    source.close()
            target.close()
    def retrieve(self):
        queue = []
        for key in self._links_dict:
            for url in self._links_dict[key]:
                queue.append({"url": url, "out_dir": os.path.join(self.reference_dir, key)})
        Utilities.single_core_queue(self.download_file, queue)
        print("Download completed")
    def merge(self):
        for fasta_type, fasta_file in zip(self._fasta_types, [self.nfasta, self.pfasta]):
            source = os.path.join(self.reference_dir, fasta_type)
            self.cat(*self.ls(source), fasta_file)
        self.index_dir = self.describer.get_index_guide(self.nfasta)
        print("Merge completed")


class Annotator:
    def __init__(self, pfasta: str):
        self.describer = ReferenceDescriber()
        self.annotation_file = self.describer.get_refdata_dict().get("sequence_1").annotation_file
        self._raw_pfasta_file = pfasta
        self._raw_nfasta_df, \
            self._processed_nfasta_df, self.nfasta_df, self.pfasta_df, self.merged_df = (pd.DataFrame(),) * 5
    @staticmethod
    def _mp_parse_nfasta_header(header):
        out = {"former_id": header}
        for tag in re.findall("\[(.+)\]", header):
            header = header.replace("[{}]".format(tag), "[{}]".format(tag.strip()))
        header_chunks = [i.strip() for i in header.split("|")]
        category_chunk = header_chunks[1].upper()
        if category_chunk.startswith("T"):
            out["category"] = "Toxin"
        elif category_chunk.startswith("AT"):
            out["category"] = "Antitoxin"
        elif category_chunk.startswith("RE"):
            out["category"] = "Regulator"
        else:
            raise ValueError("Cannot define the header's category: {}".format(header))
        out["tadb_id"] = Utilities.safe_findall("([0-9]+)", category_chunk)
        out["geninfo_id"] = Utilities.safe_findall("gi\|([0-9]+)\|", header.lower())
        ref = Utilities.safe_findall("REF\|(.+)\|", header.upper()).split("|")[0]
        if len(ref) == 0:
            try:
                ref = Utilities.safe_findall("((N|Y)(C|P|Z)_\d+\.\d*)", header.upper())[0].split("|")[0]
            except IndexError:
                pass
        out["refseq_id"] = ref
        locus = Utilities.safe_findall("\|:([c]{0,1}[0-9\-]+)", header.lower())
        out["is_antisense_strand"] = locus.startswith("c")
        out["locus"] = locus.replace("c", "")
        tail = header_chunks[-1]
        out["description"] = tail.split("[")[0].replace(
            ":{}{}".format(["", "c"][out["is_antisense_strand"]], out["locus"]), "").replace(ref, "")
        out["gene_symbol"] = Utilities.safe_findall("\[([^\[]+)\]$", tail)
        host = ""
        if tail.count("[") > 1:
            host = Utilities.safe_findall("\[([^\[]+)\]", tail, 0)
        if len(host) == 0:
            host = Utilities.safe_findall("([A-Z][a-z]+ [a-z]+[\.]{0,1})", tail)
        out["host"] = host
        for key in out:
            if isinstance(out[key], str):
                out[key] = out[key].strip()
        return out
    def annotate(self):
        # Process nucleotide FASTA
        self._raw_nfasta_df = pd.read_table(self.annotation_file, sep="\t", header=0)
        raw_nfasta_headers = self._raw_nfasta_df["former_id"].values.tolist()
        processed_nfasta_headers = [Utilities.dict2pd_series(i) for i in
                                    Utilities.multi_core_queue(self._mp_parse_nfasta_header, raw_nfasta_headers)]
        self._processed_nfasta_df = Utilities.merge_pd_series_list(processed_nfasta_headers).sort_values("former_id")
        self.nfasta_df = Utilities.left_merge(self._raw_nfasta_df, self._processed_nfasta_df, "former_id")
        # Process protein FASTA
        raw_pfasta_headers = sorted(set([j for j in [re.sub("^>", "", i).strip() for i in
                                                     open(self._raw_pfasta_file, mode="r", encoding="utf-8") if
                                                     i.startswith(">")] if len(j) > 0]))
        processed_pfasta_headers = [Utilities.dict2pd_series(i) for i in
                                    Utilities.multi_core_queue(self._mp_parse_nfasta_header, raw_pfasta_headers)]
        self.pfasta_df = Utilities.merge_pd_series_list(processed_pfasta_headers).sort_values("former_id")
        self.pfasta_df.rename(columns={"geninfo_id": "protein_geninfo_id", "refseq_id": "genpept_id",
                                       "description": "protein_description", "host": "protein_host"}, inplace=True)
        self.merged_df = Utilities.left_merge(self.nfasta_df, self.pfasta_df, "tadb_id", "category", "gene_symbol")
        self.merged_df = Utilities.combine_duplicate_rows(self.merged_df, "reference_id")
    def export(self):
        import shutil
        shutil.copy2(self.annotation_file, "{}.bak".format(self.annotation_file))
        self.merged_df.to_csv(self.annotation_file, sep='\t', index=False, header=True)


if __name__ == '__main__':
    # Paste updated version here
    retriever = SequenceRetriever("2.0")
    # retriever.retrieve()
    retriever.merge()
    """
    # Reference indexing (from worker node):
    
    rm -rf /data/reference/TADB/tadb_v2.0/index
    export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
    docker pull $IMG && \
    docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
    python3 /home/docker/scripts/cook_the_reference.py \
    -i /data/reference/TADB/tadb_v2.0/tadb_v2.0.fasta \
    -o /data/reference/TADB/tadb_v2.0/index
    
    # Wait until REFDATA file creates and complete the describer class template
    """
    retriever.describer.set_refdata("/data/reference/TADB/tadb_v2.0/index/tadb_v2.0_refdata.json")
    """
    Please update the following script lines:
    class ReferenceDescriber(ReferenceDescriberTemplate):
        NAME = "TADB"
        VERSION = "2.0"
        ALIAS = "tadb_v2.0"
        DESCRIPTION = "An updated database of bacterial type II toxin-antitoxin loci"
        DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/pubmed/?term=29106666"
        WEBSITE = "http://202.120.12.135/TADB2/index.php"
        REFDATA = "/data/reference/TADB/tadb_v2.0/index/tadb_v2.0_refdata.json"
    """
    annotator = Annotator("/data/reference/TADB/tadb_v2.0/tadb_v2.0_protein.fasta")
    annotator.annotate()
    # print(annotator.nfasta_df[annotator.nfasta_df["reference_id"].duplicated()])
    annotator.export()
