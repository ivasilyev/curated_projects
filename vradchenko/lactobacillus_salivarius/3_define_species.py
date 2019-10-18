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
import json
import pandas as pd
from copy import deepcopy
from time import sleep
from random import randint
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqRecord import SeqRecord
from meta.scripts.Utilities import Utilities
from vradchenko.lactobacillus_salivarius.ProjectDescriber import ProjectDescriber


def randomize_sleep(min_: int = 30, max_: int = 120):
    sleep(randint(min_, max_))


def attempt_func(func, args):
    _ATTEMPTS = 5
    attempt = 1
    while attempt <= _ATTEMPTS:
        try:
            if any(isinstance(args, i) for i in (list, tuple)):
                return func(*args)
            if any(isinstance(args, i) for i in (dict, )):
                return func(**args)
        except Exception as e:
            print("Caught exception for attempt {}: `{}`".format(attempt, e))
            attempt += 1
            randomize_sleep()
    print("Exceeded number of attempts for the function: '{}'".format(func.__name__))
    return


def randomize_gene_slice(record: SeqRecord):
    _LENGTH_LIMIT = 2 * (10 ** 4)
    # The typical gene is about 1000 bp in length: http://bioscience.jbpub.com/cells/MBIO137.aspx
    # The slicing will return a chunk containing ~20 genes
    gene_length = len(record)
    if gene_length <= _LENGTH_LIMIT:
        return record
    start = randint(0, gene_length - _LENGTH_LIMIT)
    end = start + _LENGTH_LIMIT
    record_ = deepcopy(record)
    record_.seq = record_.seq[start:end]
    return record_


def mp_get_and_blast_largest_contig(assembly_file: str):
    if os.path.getsize(assembly_file) == 0:
        print("Cannot process the empty file: '{}'".format(assembly_file))
        return
    with open(assembly_file) as f:
        contig_records = sorted(list(SeqIO.parse(f, "fasta")), key=lambda x: len(x), reverse=True)
        f.close()
    largest_contig = randomize_gene_slice(contig_records[0]).format("fasta")
    # The delay to avoid NCBI ban
    randomize_sleep()
    # NCBI query
    result_handle = attempt_func(NCBIWWW.qblast, ("blastn", "nt", largest_contig))
    blast_record = NCBIXML.read(result_handle)
    # Based on: https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc95
    _E_VALUE_THRESH = 0.04
    _QUERY_REPORT_SYMBOLS = 75
    high_scoring_pairs = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < _E_VALUE_THRESH:
                high_scoring_pairs.append(dict(title=alignment.title, length=alignment.length, expect=hsp.expect,
                                               score=hsp.score, bits=hsp.bits, identities=hsp.identities,
                                               positives=hsp.positives, assembly_file=assembly_file,
                                               query="...\n".join([hsp.query[:_QUERY_REPORT_SYMBOLS],
                                                                   hsp.match[:_QUERY_REPORT_SYMBOLS],
                                                                   hsp.sbjct[:_QUERY_REPORT_SYMBOLS], ""])))
    high_scoring_pairs = sorted(high_scoring_pairs, key=lambda x: x.get("score"), reverse=True)
    # Export BLAST results
    Utilities.dump_string(json.dumps(high_scoring_pairs, sort_keys=True, indent=4),
                          "{}.BLAST.json".format(os.path.splitext(assembly_file)[0]))
    return high_scoring_pairs


def process_blast_report(high_scoring_pairs: list):
    first_report = high_scoring_pairs[0]
    reference_header = first_report.get("title")
    accession_id = Utilities.safe_findall("\|* *gi\| *([^|]+) *\|", reference_header)
    return dict(assembly_file=first_report.get("assembly_file"), reference_header=reference_header,
                accession_id=accession_id)


def mp_download_reference_genbank(d: dict):
    assembly_file = d.get("assembly_file")
    accession_id = d.get("accession_id")
    Entrez.email = "name@domain.com"
    randomize_sleep()
    # NCBI query
    handle = attempt_func(Entrez.efetch, dict(db="nucleotide", id=accession_id, rettype="gb", retmode="text"))
    genbank_records = list(SeqIO.parse(handle, "genbank"))
    return dict(assembly_file=assembly_file, genbank_records=genbank_records)


def process_genbank_report(d: dict):
    genbank_records = d.get("genbank_records")
    genbank_record = genbank_records[0]
    cds_number = len([i for i in genbank_record.features if i.type == "CDS"])
    qualifiers_dict = [i.qualifiers for i in genbank_record.features if i.type == "source"][0]
    organism = Utilities.remove_empty_values(qualifiers_dict.get("organism")[0].split(" "))[:2]
    strain = " ".join(organism + [qualifiers_dict.get("strain")[0]])
    taxonomy_id = Utilities.safe_findall("\d+", [i for i in qualifiers_dict.get("db_xref") if
                                                 i.split(":")[0].strip() == "taxon"][0])
    return dict(assembly_file=d.get("assembly_file"), strain=strain, taxonomy_id=taxonomy_id,
                reference_accession_id=genbank_record.id, cds_number=cds_number, reference_bp=len(genbank_record),
                reference_description=genbank_record.description)


assemblies = [i for i in
              Utilities.scan_whole_dir(os.path.join(ProjectDescriber.ROOT_DIR, "pga-pe", "06_plasmid_merger")) if
              i.endswith(".fna") and os.path.getsize(i) > 0]
# Browse properties of the largest contigs for each assembly
props = {i: sorted(list(SeqIO.parse(i, "fasta")), key=lambda x: len(x), reverse=True)[0].format("fasta")
         for i in assemblies}
props_stats = {k: {"length": len(props.get(k)), "head": props.get(k)[:50]} for k in props}

# Create BLAST queries
blast_reports = Utilities.multi_core_queue(mp_get_and_blast_largest_contig, assemblies)
headers = Utilities.single_core_queue(process_blast_report, blast_reports)

# Create GenBank queries
genbank_reports = Utilities.multi_core_queue(mp_download_reference_genbank, headers)
reference_df = pd.DataFrame(Utilities.single_core_queue(process_genbank_report, genbank_reports))
reference_df["sample_name"] = reference_df["assembly_file"].apply(
    lambda x: "_".join(os.path.splitext(os.path.basename(x))[0].split("_")[:-1]))
reference_df.sort_values("sample_name", inplace=True)
reference_table = os.path.join(ProjectDescriber.SAMPLE_DATA_DIR, "BLASTed.sampledata")

Utilities.dump_tsv(reference_df, reference_table)
