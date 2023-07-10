
import re
from meta.utils.primitive import remove_empty_values, safe_findall


def parse_taxa(taxa):
    out = dict(genus="", species="", strain="")
    if isinstance(taxa, str):
        taxa_list = remove_empty_values(re.split("[. ]+", taxa.strip()))
        if len(taxa_list) == 0:
            return out
        out["genus"] = taxa_list[0]
        if len(taxa_list) > 1:
            if any(i.isdigit() for i in taxa_list[1]):  # Case 'Escherichia O157:H7'
                out["strain"] = taxa_list[1]
            else:  # Case 'Escherichia coli'
                sp = taxa_list[1]
                if sp in "spp":  # Case 'Escherichia sp.'
                    out["species"] = "sp."
                else:
                    out["species"] = sp
        if len(taxa_list) > 2:  # Case 'Escherichia coli O157:H7'
            out["strain"] = taxa_list[2]
    if isinstance(taxa, dict):
        out["genus"], out["species"], out["strain"] = [
            j if j is not None else ""
            for j in [
                taxa.get(i) for i in out.keys()
            ]
        ]
    return out


def regex_based_tokenization(regex_dict: dict, string: str, include_source: bool = True,
                             verbose: bool = False):
    # Column name, regex to extract, regex to excise
    out = dict()
    if include_source:
        out["source_string"] = string
    _string = str(string)

    for key, regexes in regex_dict.items():
        _string = _string.strip()
        if len(regexes) == 1:
            extract_regex = excise_regex = regexes[0]
        else:
            extract_regex = regexes[0]
            excise_regex = regexes[1]
        out[key] = safe_findall(extract_regex, _string, verbose=verbose).strip()
        _string = re.sub(excise_regex, "", _string)
    return out


def get_most_similar_pair_for_word(word: str, words: list):
    from difflib import SequenceMatcher
    # For more complicated matching, use stemming algorithms
    _words = [i for i in sorted(set(words)) if i != word]
    sorted_ratios = sorted(
        {i: SequenceMatcher(a=word, b=i).ratio() for i in _words}.items(),
        key=lambda x: x[1],
        reverse=True
    )
    try:
        out = sorted_ratios[0][0]
    except IndexError:
        print(f"Cannot process: '{word}', '{words}'")
        raise
    return [word, out]


def get_most_similar_word_pairs(words: list):
    import joblib as jb
    _words = sorted(set(words))
    out = []
    word_pairs = jb.Parallel()(
        jb.delayed(get_most_similar_pair_for_word)(i, _words) for i in _words
    )
    for word_pair in word_pairs:
        word_pair = sorted(word_pair)
        if word_pair not in out:
            out.append(word_pair)
    return out


def tokenize_reads_file_name(s: str):
    from os.path import basename
    d = regex_based_tokenization(
        {
            "extension": ["\.(.{2,8})$", "(\..{2,8})$"],  # E.g. '.fastq.gz'
            # The last segment is always '001'
            "last_segment": ["[^A-Za-z0-9]([A-Za-z0-9]+)$", "([^A-Za-z0-9][A-Za-z0-9]+)$"],
            "read_index": ["[^A-Za-z0-9](R[0-9]+)$", "([^A-Za-z0-9]R[0-9]+)$"],
            "lane_number": ["[^A-Za-z0-9](L[0-9]+)$", "([^A-Za-z0-9]L[0-9]+)$"],
            "sample_sheet_number": ["[^A-Za-z0-9](S[0-9]+)$", "([^A-Za-z0-9]S[0-9]+)$"],
            "sample_name": ["(.+)", "(.+)"],
        },
        basename(s)
    )
    d["reads_file"] = d.pop("source_string")
    if d["read_index"].endswith("1"):
        d["direction"] = "forward"
    elif d["read_index"].endswith("2"):
        d["direction"] = "reverse"
    else:
        print(f"Cannot parse read direction direction from the file '{s}'")
    return d
