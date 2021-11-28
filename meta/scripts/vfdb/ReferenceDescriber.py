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
import pandas as pd
from meta.scripts.Utilities import Utilities
from meta.scripts.reference_data import ReferenceDescriberTemplate


class ReferenceDescriber(ReferenceDescriberTemplate):
    NAME = "VFDB"
    DESCRIPTION = "A reference database for bacterial virulence factors"
    DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/pubmed/30395255"
    WEBSITE = "http://www.mgc.ac.cn/VFs/main.htm"


class SequenceRetriever:
    _DL_PAGE_URL = "http://www.mgc.ac.cn/VFs/download.htm"
    def __init__(self):
        self.describer = ReferenceDescriber()
        self.describer.VERSION = self._get_last_friday()
        self.reference_dir = os.path.join("/data/reference", self.describer.NAME, self.describer.ALIAS)
        links = [i for i in Utilities.scrap_links_from_web_page(self._DL_PAGE_URL) if i.endswith(".gz")]
        self._dl_queue = []
        for dl_link in links:
            dl_dir = self.reference_dir
            if "_pro" in dl_link:
                dl_dir = os.path.join(dl_dir, "protein")
            elif "_nt" in dl_link:
                dl_dir = os.path.join(dl_dir, "nucleotide")
            self._dl_queue.append((dl_link, dl_dir))
        self.nfasta = os.path.join(self.reference_dir, "{}.fasta".format(self.describer.ALIAS))
        self.pfasta = os.path.join(self.reference_dir, "{}_protein.fasta".format(self.describer.ALIAS))
    @staticmethod
    def _get_last_friday():
        # VFDB updates at every Friday
        import datetime
        import calendar
        last_friday = datetime.date.today()
        oneday = datetime.timedelta(days=1)
        while last_friday.weekday() != calendar.FRIDAY:
            last_friday -= oneday
        return "{:%Y.%m.%d}".format(last_friday)
    @staticmethod
    def _download_handler(input_tuple: tuple):
        os.makedirs(input_tuple[1], exist_ok=True)
        out_file = Utilities.download_file(url=input_tuple[0], out_dir=input_tuple[1])
        Utilities.decompress_file(out_file)
    def retrieve(self):
        tmp = Utilities.single_core_queue(self._download_handler, self._dl_queue)
        del tmp
        print("Download completed")
    def merge(self):
        Utilities.concatenate_files(*Utilities.scan_whole_dir(os.path.join(self.reference_dir, "nucleotide")),
                                    target_file=self.nfasta)
        Utilities.concatenate_files(*Utilities.scan_whole_dir(os.path.join(self.reference_dir, "protein")),
                                    target_file=self.pfasta)
        self.describer.get_index_guide(self.nfasta)
        print("# Merge completed. \n# Protein FASTA to annotate: '{}'\n".format(self.pfasta))


class Annotator:
    def __init__(self, pfasta: str):
        self.describer = ReferenceDescriber()
        self.annotation_file = self.describer.get_refdata_dict().get("sequence_1").annotation_file
        self._raw_pfasta_file = pfasta
        self._raw_nfasta_df = pd.DataFrame()
        self._processed_nfasta_df, self._processed_pfasta_df, self.pfasta_df, self.merged_df = (pd.DataFrame(),) * 4
    @staticmethod
    def _mp_parse_nfasta_header(header: str):
        _VFDB_REGEXES = (("vfdb_id", "^VFG(\d+)", "VFG{}"), ("gene_accession_id", "\(([^\(]+)\) ", "({}) "),
                         ("gene_symbol", "^\(([^\(]+)\) ", "({}) "), ("gene_host", "\[([^\]]+)\]$", "[{}]"),
                         ("gene_name", " \[([^\]]+)\] $", " [{}] "), ("gene_description", ".*", "{}"))
        out = {"former_id": header}
        # Spaces are important here
        for _tuple in _VFDB_REGEXES:
            key, regex, replacement = _tuple
            out[key] = Utilities.safe_findall(regex, header)
            if len(out.get(key)) > 0:
                header = header.replace(replacement.format(out.get(key)), "")
        return {k: out.get(k).strip() for k in out}
    @staticmethod
    def _mp_parse_pfasta_header(header: str):
        out = Annotator._mp_parse_nfasta_header(header)
        out = {k.replace("gene", "protein"): out.get(k) for k in out}
        out["protein_header"] = out.pop("former_id")
        return out
    def annotate(self):
        self._raw_nfasta_df = Utilities.load_tsv(self.annotation_file)
        raw_nfasta_headers = self._raw_nfasta_df["former_id"].values.tolist()
        processed_nfasta_headers = [Utilities.dict2pd_series(i) for i in
                                    Utilities.multi_core_queue(self._mp_parse_nfasta_header, raw_nfasta_headers)]
        self._processed_nfasta_df = Utilities.merge_pd_series_list(processed_nfasta_headers).sort_values("former_id")
        zf_len = len(max(self._processed_nfasta_df["vfdb_id"].values.tolist()))
        # Join table assembled from pFASTA headers
        raw_pfasta_headers = []
        with open(self._raw_pfasta_file, mode="r", encoding="utf-8") as _f:
            for _line in _f:
                if _line.startswith(">"):
                    raw_pfasta_headers.append(re.sub("^>", "", _line).strip())
            _f.close()
        raw_pfasta_headers = sorted(set([i for i in raw_pfasta_headers if len(i) > 0]))
        processed_pfasta_headers = [Utilities.dict2pd_series(i) for i in
                                    Utilities.multi_core_queue(self._mp_parse_pfasta_header, raw_pfasta_headers)]
        self._processed_pfasta_df = Utilities.merge_pd_series_list(processed_pfasta_headers).sort_values("protein_header")
        self._processed_pfasta_df["vfdb_id"] = self._processed_pfasta_df["vfdb_id"].str.zfill(zf_len)
        # Join provided table. Note the table file placed into the same dir with the merged protein FASTA file
        vfs_table_file = os.path.join(os.path.dirname(self._raw_pfasta_file), "VFs.xls")
        vfs_df = pd.read_excel(vfs_table_file, sheet_name="VFs", header=1).fillna("")
        vfs_df["vfdb_id"] = vfs_df["VFID"].str.extract("VF(\d+)")[0].str.zfill(zf_len)
        self.merged_df = pd.concat([i.set_index("vfdb_id").sort_index() for i in
                                    [self._processed_nfasta_df, self._processed_pfasta_df, vfs_df]], axis=1,
                                   sort=False).sort_index()
        self.merged_df.index.names = ["vfdb_id"]
        self.merged_df = self.merged_df.loc[self.merged_df["former_id"].str.len() > 0].reset_index()
        self.merged_df = Utilities.left_merge(self._raw_nfasta_df, self.merged_df, "former_id")
    def export(self):
        import shutil
        shutil.copy2(self.annotation_file, "{}.bak".format(self.annotation_file))
        Utilities.dump_tsv(self.merged_df, table_file=self.annotation_file)


if __name__ == '__main__':
    retriever = SequenceRetriever()
    retriever.retrieve()
    retriever.merge()
    """
    # Reference indexing (from worker node):
    
    rm -rf /data/reference/VFDB/vfdb_v2019.04.26/index
    export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
    docker pull $IMG && \
    docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it $IMG \
    python3 /home/docker/scripts/cook_the_reference.py \
    -i /data/reference/VFDB/vfdb_v2019.04.26/vfdb_v2019.04.26.fasta \
    -o /data/reference/VFDB/vfdb_v2019.04.26/index
    
    # Wait until REFDATA file creates and complete the describer class template
    
    # Merge completed
    # Protein FASTA to annotate: '/data/reference/VFDB/vfdb_v2019.04.26/vfdb_v2019.04.26_protein.fasta'
    """
    retriever.describer.set_refdata("/data/reference/VFDB/vfdb_v2019.04.26/index/vfdb_v2019.04.26_refdata.json")
    """
    Please update the following script lines:
    class ReferenceDescriber(ReferenceDescriberTemplate):
        NAME = "VFDB"
        VERSION = "2019.04.26"
        ALIAS = "vfdb_v2019.04.26"
        DESCRIPTION = "A reference database for bacterial virulence factors"
        DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/pubmed/30395255"
        WEBSITE = "http://www.mgc.ac.cn/VFs/main.htm"
        REFDATA = "/data/reference/VFDB/vfdb_v2019.04.26/index/vfdb_v2019.04.26_refdata.json"
    """
    annotator = Annotator("/data/reference/VFDB/vfdb_v2019.04.26/vfdb_v2019.04.26_protein.fasta")
    annotator.annotate()
    annotator.export()
