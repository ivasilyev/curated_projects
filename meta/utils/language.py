
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
        token = safe_findall(extract_regex, _string, verbose=verbose)


        out[key] = safe_findall(extract_regex, _string, verbose=verbose).strip()
        _string = re.sub(excise_regex, "", _string)
    return out
