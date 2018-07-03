#!/usr/bin/env python3
# -*- coding: utf-8 -*-



# Ex.: https://www.ncbi.nlm.nih.gov/protein/?term="Escherichia+coli"[orgn]+AND+"alpha+hemolysin"


import re
import requests
import bs4
import pandas as pd
from meta.scripts.utilities import dict2pd_series
from meta.scripts.FASTAParser import FASTA


class ProteinSequenceRetriever:
    """
    This class performs NCBI Protein DB search. Consumes organism name (space-delimited) and gene name.
    Query example: 'https://www.ncbi.nlm.nih.gov/protein/?term="Escherichia+coli"[orgn]+AND+"alpha+hemolysin"'
    """
