
def get_headers_from_fasta(file: str):
    from re import sub
    from meta.utils.primitive import remove_empty_values
    out = []
    with open(file, mode="r", encoding="utf-8") as _f:
        for _line in _f:
            if _line.startswith(">"):
                out.append(sub("^>", "", _line).strip())
        _f.close()
    return sorted(set(remove_empty_values(out)))
