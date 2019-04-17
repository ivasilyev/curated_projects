#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export DOCKER_IMAGE_NAME=ivasilyev/curated_projects:latest && \
docker pull ${DOCKER_IMAGE_NAME} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${DOCKER_IMAGE_NAME} python3
"""

import re
import pandas as pd
import shutil
from meta.templates.ReferenceDescriberTemplate import ReferenceDescriberTemplate
from meta.scripts.Utilities import Utilities


class ReferenceDescriber(ReferenceDescriberTemplate):
    alias = "tadb_v2.0"
    name = "TADB"
    description = "An updated database of bacterial type II toxin-antitoxin loci"
    documentation = "https://www.ncbi.nlm.nih.gov/pubmed/?term=29106666"
    url = "http://202.120.12.135/TADB2/index.php"
    refdata = "/data/reference/TADB/index/tadb_v2.0/tadb_v2.0_refdata.json"


class Annotator:
    def __init__(self):
        self._describer = ReferenceDescriber()
        self._annotation_file = self._describer.get_refdata_dict().get("sequence_1").annotation_file
        self.nfasta_df = None
        self.pfasta_df = None
    def _mp_parse_nfasta_header(self, header):
        output_dict = {"former_id": header}
        for tag in re.findall("\[(.+)\]", header):
            header = header.replace("[{}]".format(tag), "[{}]".format(tag.strip()))
        header_chunks = [i.strip() for i in header.split("|")]
        category_chunk = header_chunks[1].upper()
        if category_chunk.startswith("T"):
            output_dict["category"] = "Toxin"
        elif category_chunk.startswith("AT"):
            output_dict["category"] = "Antitoxin"
        elif category_chunk.startswith("RE"):
            output_dict["category"] = "Regulator"
        else:
            raise ValueError("Cannot define the header's category: {}".format(header))
        output_dict["tadb_id"] = Utilities.safe_findall("([0-9]+)", category_chunk)
        output_dict["geninfo_id"] = Utilities.safe_findall("gi\|([0-9]+)\|", header.lower())
        output_dict["genbank_id"] = Utilities.safe_findall("REF\|(.+)\|", header.upper()).split("|")[0]
        output_dict["locus"] = Utilities.safe_findall("\|:([c]{0,1}[0-9\-]+)", header.lower())
        output_dict["is_antisense_strand"] = output_dict["locus"].startswith("c")
        output_dict["host"] = Utilities.safe_findall("([A-Z][a-z]* [a-z]+[\.]{0,1})", header_chunks[-1])
        output_dict["gene_name"] = Utilities.safe_findall("\[(.+)\]$", header_chunks[-1])
        return Utilities.dict2pd_series(output_dict)
    def annotate(self):
        self.nfasta_df = pd.read_table(self._annotation_file, sep="\t", header=0).set_index("former_id")
        mp_queue = self.nfasta_df.index.values.tolist()
        mp_result = Utilities.multi_core_queue(self._mp_parse_nfasta_header, mp_queue)
        mp_df = Utilities.merge_pd_series_list(mp_result).set_index("former_id").sort_index()
        shutil.copy2(self._annotation_file, "{}.bak".format(self._annotation_file))
        self.nfasta_df = pd.concat([self.nfasta_df, mp_df], axis=1, sort=False)


if __name__ == '__main__':
    # TODO if db becomes rapidly updating: auto crawl download page: http://202.120.12.135/TADB2/download.html
    dl_list = [
        # The experimentally validated Type II TA loci
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_T_exp.fas",  # Type II Toxin
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_AT_exp.fas",  # Type II Antitoxin (protein)
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_RE.fas",  # Type II Regulator (of 3-component type II TA loci)
        # All the in silico predicted and experimentally validated Type II TA loci
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_T.fas",
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_AT.fas",
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_II_RE.fas",
        # The experimentally validated Type I TA loci
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_I_T_exp.fas",
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_I_AT_exp.fas",
        # All the in silico predicted and experimentally validated Type I TA loci
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_I_T.fas",
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_I_AT.fas",
        # The experimentally validated Type III TA loci
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_III_T.fas",
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_III_AT.fas",
        # The experimentally validated Type IV TA loci
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_IV_T.fas",
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_IV_AT.fas",
        # The experimentally validated Type V TA loci
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_V_T.fas",
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_V_AT.fas",
        # The experimentally validated Type VI TA loci
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_VI_T.fas",
        "http://202.120.12.135/TADB2/download/TADB2/20171013/nucleotide/type_VI_AT.fas"
    ]

    describer = ReferenceDescriber()
    describer.get_index_guide("/data/reference/TADB/tadb_v2.0.fasta")
    annotator = Annotator()
    annotator.annotate()

"""
# Reference indexing (from worker node):

rm -rf /data/reference/TADB/index/tadb_v2.0
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull $IMG && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
python3 /home/docker/scripts/cook_the_reference.py \
-i /data/reference/TADB/tadb_v2.0.fasta \
-o /data/reference/TADB/index/tadb_v2.0

# Wait until REFDATA file creates and complete the describer class template
"""
