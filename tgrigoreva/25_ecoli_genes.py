#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests
import bs4
import re
import multiprocessing

# Pre-setup:
# docker pull debian:jessie && docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it debian:jessie bash
# apt-get -y update && apt-get -y install ca-certificates curl python3 python3-pip


class NuccoreSequenceRetriever:
    """
    This class performs NCBI Gene DB search. Consumes organism name (space-delimited) and gene name.
    May return raw pages or whole FASTA from valid pages
    """
    def __init__(self, species, gene):
        self._species = species.strip()
        self._gene = gene.strip()
        self._headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.9; rv:45.0) Gecko/20100101 Firefox/45.0'}
        self._search_url = "https://www.ncbi.nlm.nih.gov/gene/?term=%22" + self._species.replace(' ', '+') + "%22%5Borgn%5D+AND+" + self._gene
        self._search_soup = bs4.BeautifulSoup(requests.get(self._search_url, headers=self._headers).content, "lxml")
        self._row_soups_list = self._search_soup.find_all("tr", "rprt")
    def get_soup(self):
        return self._search_soup
    @staticmethod
    def parse_table_row(row_soup):
        d = {"Name": row_soup.find_all("td", "gene-name-id")[0].find_all("a")[0].text,
             "Gene ID": "".join(re.findall("ID: (\d+)",row_soup.find_all("td", "gene-name-id")[0].find_all("span", "gene-id")[0].text)),
             "Description": row_soup.find_all("td")[1].text,
             "Sequence ID": "".join(re.findall("(NC_\d+\.\d*)", row_soup.find_all("td")[2].text)),
             "Sequence Location": "".join(re.findall("\((\d+\.\.\d+)", row_soup.find_all("td")[2].text)),
             "Aliases": row_soup.find_all("td")[3].text}
        return {k: d[k].strip() for k in d}
    def _get_fasta(self, organism_id, coordinates_list):
        if not isinstance(coordinates_list, list):
            raise ValueError("Coordinates must be in a list")
        coordinates_list = [str(i) for i in coordinates_list]
        _soup = bs4.BeautifulSoup(requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + organism_id + "&rettype=fasta&retmode=text&seq_start=" + coordinates_list[0] + "&seq_stop=" + coordinates_list[-1], headers=self._headers).content, "lxml")
        _fasta = re.sub("\n+", "\n", _soup.find_all("p")[0].text)
        return _fasta
    def query2fasta(self):
        def _process_fasta_header(_fasta):
            _fasta_sequence_list = _fasta.split('\n')[1:]
            _fasta_header_string = _fasta.split('\n')[0] + " RETRIEVED BY " + self._gene + " WITH NANE " + _row_dict["Name"]
            return "\n".join([_fasta_header_string] + _fasta_sequence_list)
        out = []
        for _soup in self._row_soups_list:
            _row_dict = self.parse_table_row(_soup)
            _locations_list = _row_dict["Sequence Location"].split("..")
            # Filtering expression
            # if len(_locations_list) == 2:
            if (self._gene.lower() in _row_dict["Name"].lower() or self._gene.lower() in _row_dict["Description"].lower() or self._gene.lower() in _row_dict["Aliases"].lower()) and len(_locations_list) == 2:
                fasta = self._get_fasta(_row_dict["Sequence ID"], _locations_list)
                out.append(_process_fasta_header(fasta))
        return "\n".join(out)


genesList = "Stx2, EhxA, STb, EspA, EspB, EspC, Cnf, Cfa, Iha, pap, papA, papC, papE, papF, Tir, Etp, KpsM, KpsT, FliC, IbeA, Tsh, IucD, TraT, IutA, espP, katP, ompA, ompT, iroN, iss, fyuA, uidA, uspA, cdtB, cvaC, ibeA".split(", ")


def multi_core_queue(function_to_parallelize, queue):
    pool = multiprocessing.Pool()
    output = pool.map(function_to_parallelize, queue)
    pool.close()
    pool.join()
    return output


def get_fasta_by_gene(name):
    obj = NuccoreSequenceRetriever("Escherichia coli", name)
    return obj.query2fasta()


def process_genes_list():
    return multi_core_queue(get_fasta_by_gene, genesList)


genesSequencesQueueList = process_genes_list()


def filter_sequences_list():
    d = {}
    for i in genesSequencesQueueList:
        for j in i.split("\n\n"):
            if len(j) > 0:
                k = re.findall(">(.*) RETRIEVED", j)[0]
                if k not in d:
                    d.update({k: j})
    return [d[k].strip() for k in d]


genesSequencesFilteredList = filter_sequences_list()
