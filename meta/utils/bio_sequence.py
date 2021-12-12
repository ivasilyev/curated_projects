
def remove_duplicate_sequences(records: list):
    out = []
    sequences = []
    for record in records:
        if record.seq not in sequences:
            sequences.append(record.seq)
            out.append(record)
    return out


def load_sequences(file: str, fmt: str = "fasta"):
    from Bio import SeqIO
    from meta.utils.primitive import remove_empty_values
    if fmt == "fastq_gz":
        import gzip
        with gzip.open(file, "rt") as f:
            records = list(SeqIO.parse(f, "fastq"))
            f.close()
    else:
        with open(file, mode="r", encoding="utf-8") as f:
            records = list(SeqIO.parse(f, fmt))
            f.close()
    out = remove_empty_values(sorted(records, key=lambda x: len(x), reverse=True))
    return out


def string_to_sequences(s: str, fmt: str = "fasta"):
    from Bio import SeqIO
    from io import StringIO
    return list(SeqIO.parse(StringIO(s), fmt))


def dump_sequences(sequences: list, file: str, fmt: str = "fasta"):
    import os
    from Bio import SeqIO
    os.makedirs(os.path.dirname(file), exist_ok=True)
    with open(file, mode="w", encoding="utf-8") as f:
        SeqIO.write(sequences, f, fmt)
        f.close()


def load_headers_from_fasta(file: str):
    from re import sub
    from meta.utils.primitive import remove_empty_values
    out = []
    with open(file, mode="r", encoding="utf-8") as _f:
        for _line in _f:
            if _line.startswith(">"):
                out.append(sub("^>", "", _line).strip())
        _f.close()
    return sorted(set(remove_empty_values(out)))


def randomize_gene_slice(record, size: int = 20000):
    """
    :param record: SeqRecord

    :param size: int

    Cuts a slice with random start and of given length from sequence record
    The typical gene is about 1000 bp in length:
    http://bioscience.jbpub.com/cells/MBIO137.aspx
    The default slicing will return a chunk containing ~20 genes.

    :return: trimmed SeqRecord
    """
    from copy import deepcopy
    from random import randint
    gene_length = len(record)
    if gene_length <= size:
        return record
    start = randint(0, gene_length - size)
    end = start + size
    record_ = deepcopy(record)
    record_.seq = record_.seq[start:end]
    return record_


def describe_genbank(genbank_record, verbose: bool = True):
    from Bio.SeqUtils import GC
    from Bio.GenBank import Record
    assert type(genbank_record) == Record
    out = dict(
        total_cds=0,
        genbank_id=genbank_record.id,
        gc_percentage=round(GC(genbank_record.seq), 2),
        reference_bp=len(genbank_record),
        reference_description=genbank_record.description
    )
    try:
        out["total_cds"] = int(
            genbank_record.annotations["structured_comment"]["Genome-Annotation-Data"]["CDSs (total)"].replace(",", "")
        )
    except KeyError:
        pass
    if out["total_cds"] == 0:
        out["total_cds"] = len([i for i in genbank_record.features if i.type == "CDS"])
    if out["total_cds"] == 0 and verbose:
        print("Cannot find annotations for GenBank Accession ID '{}'".format(out["genbank_id"]))
    annotation_dict = {k: v for k, v in genbank_record.annotations.items()
                       if all(not isinstance(v, i) for i in [dict, list])}
    out.update(annotation_dict)
    if "taxonomy" in genbank_record.annotations.keys():
        out["taxonomy"] = ";".join(genbank_record.annotations["taxonomy"])
    qualifiers_dict = [i.qualifiers for i in genbank_record.features if i.type == "source"][0]
    out.update({k: ";".join(v) for k, v in qualifiers_dict.items()})
    if "db_xref" in out.keys():
        out["taxonomy_id"] = out["db_xref"].replace("taxon:", "")
    return out

