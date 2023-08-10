
import os
import pandas as pd
from meta.utils.io import load_dict
from meta.utils.web import get_file
from meta.utils.pandas import dump_tsv
from meta.utils.date_time import get_timestamp

reference_db_name = "KEGG"
output_dir = "/data/reference"

timestamp = get_timestamp(fmt="%Y-%m-%d")
reference_mask = f"{reference_db_name.lower()}_v{timestamp}"
reference_dir = os.path.join(output_dir, reference_db_name, reference_mask)

reference_prefix = os.path.join(reference_dir, f"{reference_mask}")
reference_json = f"{reference_prefix}.json"

get_file(
    "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir=",
    reference_json,
    force=False
)
reference_dict = load_dict(reference_json)

out_dicts = list()


def _expand(d: dict, names: list = None):
    if not isinstance(names, list):
        names = []
    name = d["name"]
    _names = names
    if name != "ko00001":
        _names = names + [name]
    if "children" in d.keys():
        for child in d["children"]:
            _expand(child, _names)
    else:
        _name = name.split(" ")[0]
        out_dicts.append(dict(ko=_name, **{f"L{idx + 1}": i for idx, i in enumerate(_names)}))


_expand(reference_dict)

reference_df = pd.DataFrame(out_dicts).set_index("ko").sort_index()
dump_tsv(reference_df, f"{reference_prefix}_denormalized.tsv", reset_index=True)
