#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import re
import requests
import bs4
import pandas as pd
from meta.scripts.utilities import dict2pd_series
from meta.scripts.FASTAParser import FASTA


class NuccoreSequenceRetriever:
    """
    This class performs NCBI Gene DB search. Consumes organism name (space-delimited) and gene name.
    Query example: 'https://www.ncbi.nlm.nih.gov/gene/?term="Escherichia+coli"[orgn]+AND+"alpha+hemolysin"'
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
    def parse_table_row(self, row_soup):
        d = {"Name": row_soup.find_all("td", "gene-name-id")[0].find_all("a")[0].text,
             "Gene ID": "".join(re.findall("ID: (\d+)", row_soup.find_all("td", "gene-name-id")[0].find_all("span", "gene-id")[0].text)),
             "Description": row_soup.find_all("td")[1].text,
             "Sequence ID": "".join(re.findall("(NC_\d+\.\d*)", row_soup.find_all("td")[2].text)),
             "Sequence Location": "".join(re.findall("\((\d+\.\.\d+)", row_soup.find_all("td")[2].text)),
             "Aliases": row_soup.find_all("td")[3].text,
             "Query": self._gene}
        return {k: d[k].strip() for k in d}
    def _get_fasta(self, organism_id, coordinates_list):
        if not isinstance(coordinates_list, list):
            raise ValueError("Coordinates must be in a list")
        coordinates_list = [str(i) for i in coordinates_list]
        _soup = bs4.BeautifulSoup(requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + organism_id + "&rettype=fasta&retmode=text&seq_start=" + coordinates_list[0] + "&seq_stop=" + coordinates_list[-1], headers=self._headers).content, "lxml")
        _fasta = FASTA(re.sub("\n+", "\n", _soup.find_all("p")[0].text))
        return _fasta
    def query2fasta(self):
        output_dict = {"FASTAs_list": [], "annotations_series_list": []}
        for _soup in self._row_soups_list:
            _row_dict = self.parse_table_row(_soup)
            _locations_list = _row_dict["Sequence Location"].split("..")
            # Filtering expression
            if (self._gene.lower() in _row_dict["Name"].lower() or self._gene.lower() in _row_dict["Description"].lower() or self._gene.lower() in _row_dict["Aliases"].lower()) and len(_locations_list) == 2:
                fasta = self._get_fasta(_row_dict["Sequence ID"], _locations_list)
                output_dict["FASTAs_list"].append(fasta)
                output_dict["annotations_series_list"].append(pd.Series(dict2pd_series(_row_dict), name=fasta.header))
        return output_dict
