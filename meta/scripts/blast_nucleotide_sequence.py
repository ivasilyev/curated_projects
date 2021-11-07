#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from Bio.SeqUtils import GC
from Bio import SeqIO, Entrez
from time import perf_counter
from collections import OrderedDict
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.GenBank import Record as GBRecord
from meta.scripts.Utilities import Utilities


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
    return (_namespace.input,
            _namespace.blast_only,
            _namespace.chromosomes_only,
            _namespace.results,
            _namespace.sequence_dir,
            _namespace.output)


def parse_largest_subsequence(fasta_nt_file: str):
    """
    Parses nucleotide FASTA and chunk if it's too large for BLAST query
    """
    assert Utilities.is_file_valid(fasta_nt_file)
    records = Utilities.parse_sequences(fasta_nt_file)
    return Utilities.randomize_gene_slice(records[0]).format("fasta")


def download_nt_blast_report(query: str):
    # The delay to avoid NCBI ban
    Utilities.randomize_sleep()
    # NCBI query
    result_handle = Utilities.attempt_func(NCBIWWW.qblast, ("blastn", "nt", query))
    return NCBIXML.read(result_handle)


def parse_blast_report(blast_record):
    _E_VALUE_THRESH = 0.04
    # Based on: https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc95
    high_scoring_pairs = OrderedDict()
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < _E_VALUE_THRESH:
                t = alignment.title
                d = dict(title=t, length=alignment.length, expect=hsp.expect, score=hsp.score,
                         bits=hsp.bits, identities=hsp.identities, positives=hsp.positives,
                         query="\n".join([hsp.query, hsp.match, hsp.sbjct, ""]),
                         geninfo_id=Utilities.safe_findall("\|* *gi\| *([^|]+) *\|", t))
                high_scoring_pairs[t] = d
    return high_scoring_pairs


def download_reference_genbank(accession_id: str):
    Entrez.email = "name@domain.com"
    # The delay to avoid NCBI ban
    Utilities.randomize_sleep(10, 30)
    # NCBI query
    handle = Utilities.attempt_func(Entrez.efetch, dict(db="nucleotide", id=accession_id,
                                                        rettype="gb", retmode="text"))
    return list(SeqIO.parse(handle, "genbank"))


def describe_reference_genbank(genbank_record: GBRecord):
    cds_number = 0
    try:
        cds_number = int(genbank_record.annotations["structured_comment"]["Genome-Annotation-Data"]["CDSs (total)"].replace(",", ""))
    except KeyError:
        pass
    if cds_number == 0:
        cds_number = len([i for i in genbank_record.features if i.type == "CDS"])
        # Otherwise the reference genbank might not be annotated
    qualifiers_dict = [i.qualifiers for i in genbank_record.features if i.type == "source"][0]
    organism = Utilities.remove_empty_values(qualifiers_dict.get("organism")[0].split(" "))[:2]
    strain = " ".join(organism + [qualifiers_dict.get("strain")[0]])
    taxonomy_id = Utilities.safe_findall("\d+", [i for i in qualifiers_dict.get("db_xref") if
                                                 i.split(":")[0].strip() == "taxon"][0])
    gc_percentage = round(GC(genbank_record.seq), 2)
    return dict(strain=strain, taxonomy_id=taxonomy_id, genbank_id=genbank_record.id,
                total_cds=cds_number, reference_bp=len(genbank_record),
                reference_description=genbank_record.description, gc_percentage=gc_percentage)


if __name__ == '__main__':
    nt_fasta_file, is_blast_only, is_chromosomes_only, blast_result_number, sequence_directory, output_directory = _parse_args()
    out_blast_basename = os.path.join(output_directory, Utilities.filename_only(nt_fasta_file))
    blast_result_file = "{}_blast_results.json".format(out_blast_basename)

    if Utilities.is_file_valid(blast_result_file, report=False):
        print("Load results of the already performed NCBI BLAST query")
        blast_results = Utilities.load_dict(blast_result_file)
    else:
        print("Parsing largest subsequence")
        blast_query = parse_largest_subsequence(nt_fasta_file)
        Utilities.dump_string(blast_query, "{}_blast_query.fna".format(out_blast_basename))
        print(f"Performing BLAST query from the sequence of length {len(blast_query)}")
        start = perf_counter()
        blast_report = download_nt_blast_report(blast_query)
        print(f"The NCBI query was completed with {Utilities.count_elapsed_seconds(start)}")
        blast_results = parse_blast_report(blast_report)
        Utilities.dump_dict(blast_results, blast_result_file, sort_keys=False, indent=4)

    if not is_blast_only:
        if len(sequence_directory) == 0:
            sequence_directory = os.path.join(output_directory, "genbank")

        genbank_descriptions = []
        counter = 1
        for blast_result in list(blast_results.values()):
            geninfo_accession = blast_result["geninfo_id"]
            genbank_file = os.path.join(sequence_directory, "{}.gbk".format(geninfo_accession))
            if Utilities.is_file_valid(genbank_file, report=False):
                genbank_report = Utilities.parse_sequences(genbank_file, "genbank")
            else:
                if is_chromosomes_only and " chromosome" not in blast_result["title"]:
                    continue
                genbank_report = download_reference_genbank(geninfo_accession)
                Utilities.dump_string(genbank_report[0].format("genbank"), genbank_file)
            genbank_description = describe_reference_genbank(genbank_report[0])
            genbank_description["geninfo_id"] = geninfo_accession
            genbank_description["genbank_file"] = genbank_file
            genbank_descriptions.append(genbank_description)
            print(f"Downloaded {counter} of {blast_result_number} sequences")
            counter += 1
            if counter == blast_result_number:
                break

        blast_result_df = pd.DataFrame(blast_results.values())
        genbank_description_df = pd.DataFrame(genbank_descriptions)
        combined_blast_result_df = pd.concat([i.set_index("geninfo_id") for i in (blast_result_df, genbank_description_df)],
                                             axis=1, sort=False).drop(["query", "genbank_file"], axis=1).dropna(how="any").rename_axis(index="geninfo_id")
        Utilities.dump_tsv(combined_blast_result_df.reset_index(), os.path.join(output_directory, "combined_blast_results.tsv"))

        report_dict = dict(input_file=nt_fasta_file, genbank_files=genbank_description_df["genbank_file"].values.tolist())
        Utilities.dump_dict(report_dict, os.path.join(output_directory, "report.json"))
    print("Completed")
