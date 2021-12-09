
import re
from meta.utils.primitive import safe_findall


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
