#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Pre-setup:
export IMG=ivasilyev/curated_projects:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 --net=host -it ${IMG} bash

git pull
python3
"""

import os
import re
import pandas as pd
from Bio import Phylo, SeqIO
from meta.scripts.Utilities import Utilities


def flatten_string(s: str):
    if not s:
        return ""
    return re.sub(" +", "",  re.sub("[\r\n ]+", " ", s))


def remove_maintenance_comments(s: str):
    _RGX = ("REFSEQ INFORMATION:", "The reference sequence was derived from [^\.]+\.", "https:[^ ]+",
            "The annotation was added by the NCBI Prokaryotic Genome Annotation Pipeline \(PGAP\)\.",
            "Information about PGAP can be found here:",
            "Annotation was added by the NCBI Prokaryotic Genome Annotation Pipeline \(released 2013\)\.",
            "Information about the Pipeline can be found here:",
            "Annotation Pipeline \(PGAP\) set[;]{0,1}", "Annotation Pipeline set[;]{0,1}",
            "GeneMark[^ ]+ repeat_region[;]{0,1}", "COMPLETENESS: full length\.")
    for regex in _RGX:
        s = re.sub(regex, "", s)
    return flatten_string(s)


genbank_dir = "/data1/bio/projects/inicolaeva/klebsiella_infants/ncbi-dl/gbff"
genbank_files = [i for i in Utilities.scan_whole_dir(genbank_dir) if i.endswith(".gbff")]

tree = Phylo.read("/data1/bio/projects/inicolaeva/klebsiella_infants/roary/newick/iTOL_collapsed_tree.newick", "newick")
node_names = [j for j in [i.name for i in tree.find_clades()] if j is not None and j.startswith("GCF")]

annotations_list = []
for node_name in node_names:
    # node_name = "GCF_005377825.1_ASM537782v1"
    genbank_file = os.path.join(genbank_dir, "{}_genomic.gbff".format(node_name))
    seq_records = list(SeqIO.parse(genbank_file, "genbank"))
    annotation_dict = {i: flatten_string(seq_records[0].annotations.get(i)) for i in ["organism", "date", "comment"]}
    annotation_dict["comment"] = remove_maintenance_comments(annotation_dict["comment"])
    annotation_dict["strain"] = Utilities.safe_findall("[S|s]train:* ([^ ]+)", seq_records[0].description)
    annotation_dict["refseq_id"] = Utilities.safe_findall("GCF_[^_]+", node_name)
    annotation_dict["assembly_id"] = node_name.replace(annotation_dict["refseq_id"], "").strip("_")
    annotations_list.append(annotation_dict)

annotations_df = pd.DataFrame(annotations_list)
Utilities.dump_tsv(annotations_df,
                   "/data1/bio/projects/inicolaeva/klebsiella_infants/roary/newick/iTOL_collapsed_tree_annotation.tsv")
