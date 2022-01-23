#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from time import perf_counter
from collections import OrderedDict
from Bio import SeqIO, SeqRecord, Entrez
from Bio.Blast import NCBIWWW, NCBIXML, Record
from meta.utils.queue import attempt_func
from meta.utils.pandas import dump_tsv, concat
from meta.utils.io import load_dict, dump_dict, dump_string
from meta.utils.file_system import is_file_valid, filename_only
from meta.utils.primitive import safe_findall, clear_non_printing_chars
from meta.utils.date_time import randomize_sleep, count_elapsed_seconds
from meta.utils.bio_sequence import describe_genbank, dump_sequences, load_sequences, randomize_gene_slice


E_VALUE_THRESH = 0.04
SLEEP_INTERVAL = (5, 15)
QUERY_SIZE = 20000
QUERY_ATTEMPTS = 5


def download_nt_blast_report(query: str, result_number: int = 50):
    # The delay to avoid NCBI ban
    randomize_sleep(*SLEEP_INTERVAL)
    # NCBI query
    result_handle = attempt_func(
        NCBIWWW.qblast, database="nt", program="blastn", hitlist_size=result_number, sequence=query
    )
    return NCBIXML.read(result_handle)


def parse_blast_report(blast_record: Record):
    # Based on: https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc95
    high_scoring_pairs = OrderedDict()
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                t = alignment.title
                d = dict(
                    title=t,
                    length=alignment.length,
                    expect=hsp.expect,
                    score=hsp.score,
                    bits=hsp.bits,
                    identities=hsp.identities,
                    positives=hsp.positives,
                    query=hsp.query,
                    match=hsp.match,
                    sbjct=hsp.sbjct,
                    geninfo_id=safe_findall("gi\|([^|]+)\|", t).strip()
                )
                high_scoring_pairs[t] = {k: clear_non_printing_chars(v) for k, v in d.items()}
    return high_scoring_pairs


def chop_and_blast(record: SeqRecord, chunk_size: int = QUERY_SIZE, result_number: int = 50):
    print(f"Getting random subsequence chunk of size {chunk_size} from the sequnce of size {len(record)} bp.")
    query_string = randomize_gene_slice(record, size=chunk_size).format("fasta")
    print(f"Performing BLAST query from the sequence of length {len(query_string)}")
    _start = perf_counter()
    report = download_nt_blast_report(query_string, result_number)
    print(f"BLAST query was completed after {count_elapsed_seconds(_start)}")
    return query_string, parse_blast_report(report)


def download_reference_genbank(accession_id: str):
    Entrez.email = "name@domain.com"
    # The delay to avoid NCBI ban
    randomize_sleep(*SLEEP_INTERVAL)
    # NCBI query
    handle = attempt_func(
        Entrez.efetch, db="nucleotide", id=accession_id, rettype="gb", retmode="text"
    )
    return list(SeqIO.parse(handle, "genbank"))


def _parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description="Tool to perform BLAST using randomized slice from the largest nucleotide "
                    "sequence and download the related GenBank reference entries",
        epilog="Sequences must be in FASTA format")
    parser.add_argument("-i", "--input", metavar="<file>", required=True, help="FASTA nucleotide file")
    parser.add_argument("-b", "--blast_only", default=False, action="store_true",
                        help="(Optional) If selected, the related GenBank reference entries won't be downloaded")
    parser.add_argument("-c", "--chromosomes_only", default=False, action="store_true",
                        help="(Optional) If selected, only chromosome GenBank reference entries will be downloaded")
    parser.add_argument("-r", "--results", metavar="<int>", type=int, default=25,
                        help="(Optional) The maximum size of GenBank entries from the BLAST report to download")
    parser.add_argument("-s", "--sequence_dir", metavar="<directory>", default="",
                        help="(Optional) Special dir to download sequences to ('genbank' subdirectory by default)")
    parser.add_argument("-o", "--output", metavar="<directory>", required=True, help="Output directory")
    _namespace = parser.parse_args()
    return (
        _namespace.input,
        _namespace.blast_only,
        _namespace.chromosomes_only,
        _namespace.results,
        _namespace.sequence_dir,
        _namespace.output
    )


if __name__ == '__main__':
    (
        nt_fasta_file,
        is_blast_only,
        is_chromosomes_only,
        blast_result_number,
        sequence_directory,
        output_directory
    ) = _parse_args()
    out_blast_basename = os.path.join(output_directory, filename_only(nt_fasta_file))
    blast_query_file = f"{out_blast_basename}_blast_query.fna"
    blast_result_file = "{}_blast_results.json".format(out_blast_basename)

    blast_results = dict()

    if is_file_valid(blast_result_file):
        print(f"Load results of the already performed NCBI BLAST query: '{blast_result_file}'")
        blast_results = load_dict(blast_result_file)
    else:
        fasta_record = load_sequences(nt_fasta_file, "fasta")[0]
        for attempt in range(QUERY_ATTEMPTS):
            blast_query_string, blast_results = chop_and_blast(
                record=fasta_record, result_number=blast_result_number
            )
            if len(blast_results.keys()) > 0:
                dump_string(blast_query_string, blast_query_file)
                print(f"Saved BLAST query to file '{blast_query_file}'")
                dump_dict(blast_results, blast_result_file)
                print(f"{len(blast_results.keys())} BLAST results were saved into {blast_result_file}")
                break
            print(f"Empty BLAST results for attempt {attempt} of {QUERY_ATTEMPTS}")

    genbank_description_file = "{}_genbank_descriptions.json".format(out_blast_basename)
    if not is_blast_only:
        if len(sequence_directory) == 0:
            sequence_directory = os.path.join(output_directory, "genbank")

        genbank_descriptions = dict()
        for counter, (blast_result_title, blast_result) in enumerate(list(blast_results.items())):
            geninfo_accession = blast_result["geninfo_id"]
            genbank_file = os.path.join(sequence_directory, "{}.gbk".format(geninfo_accession))
            if is_file_valid(genbank_file):
                genbank_reports = load_sequences(genbank_file, "genbank")
            else:
                if is_chromosomes_only and all(
                    f" {i}" not in blast_result_title for i in (
                        "chromosome", "complete", "genome"
                    )
                ):
                    continue
                genbank_reports = download_reference_genbank(geninfo_accession)
                dump_sequences(genbank_reports, genbank_file, "genbank")

            genbank_description = describe_genbank(genbank_reports[0])
            genbank_description.update(dict(
                geninfo_id=geninfo_accession,
                genbank_file=genbank_file
            ))
            genbank_descriptions[blast_result_title] = genbank_description

        dump_dict(genbank_descriptions, genbank_description_file)
        print(f"{len(genbank_descriptions.keys())} GenBank descriptions were saved into {genbank_description_file}")

        combined_blast_result_df = concat(
            [pd.DataFrame(i.values()) for i in (blast_results, genbank_descriptions)],
            axis=1,
            how="outer",
            on="geninfo_id",
        )
        print(f"Merged {len(blast_results.keys())} BLAST results and {len(genbank_descriptions.keys())} result descriptions")
        report_dict = dict(
            input_file=nt_fasta_file,
            genbank_files=combined_blast_result_df.loc[
                combined_blast_result_df["genbank_file"].notnull() &
                combined_blast_result_df["taxonomy"].notnull(),
                :
            ]["genbank_file"].values.tolist()
        )
        report_json = os.path.join(output_directory, "report.json")
        dump_dict(report_dict, report_json)
        print("Saved report JSON of {} items into '{}'".format(len(report_dict["genbank_files"]), report_json))

        out_blast_result_df = combined_blast_result_df.drop(
            ["query", "genbank_file", "match", "sbjct"],
            axis=1
        )
        combined_blast_result_file = os.path.join(output_directory, "combined_blast_results.tsv")
        dump_tsv(out_blast_result_df, combined_blast_result_file)
        print(f"Saved combined BLAST result table with shape of {out_blast_result_df.shape} into '{combined_blast_result_file}'")
    print("Completed")
