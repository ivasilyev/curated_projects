#%%

import os
import re
import pandas as pd
from json import loads
from meta.utils.date_time import get_timestamp
from meta.utils.diversity import TAXONOMY_CHARS
from meta.utils.io import dump_dict
from meta.utils.pandas import dump_tsv
from meta.utils.primitive import remove_empty_values
from meta.utils.queue import multi_core_queue2
from meta.utils.web import get_page, get_soup, parse_links_from_soup


KO_REGEX = "ko[0-9]+$"
REFERENCE_DB_NAME = "KEGG"


def _expand_kegg_dict(d: dict, names: list = None):
    # Recursive depth-first search
    if not isinstance(names, list):
        names = []
    name = d["name"]
    _names = names
    if re.match(KO_REGEX, name) is None:  # The root node
        _names = names + [name]
    if "children" in d.keys():
        for child in d["children"]:
            _expand_kegg_dict(child, _names)
    else:
        _name = name.split(" ")[0]
        out_dicts.append(dict(ko=_name, **{f"L{idx + 1}": i for idx, i in enumerate(_names)}))


def _get_dict(ko: str):
    url = f"https://www.kegg.jp/kegg-bin/download_htext?htext={ko}.keg&format=json&filedir="
    page = get_page(url)
    s = page.decode("utf-8")
    return loads(s)

out_dicts = list()

#%%


output_dir = "/data/reference"

timestamp = get_timestamp(fmt="%Y-%m-%d")
reference_mask = f"{REFERENCE_DB_NAME.lower()}_v{timestamp}"
reference_dir = os.path.join(output_dir, REFERENCE_DB_NAME, reference_mask)

reference_prefix = os.path.join(reference_dir, f"{reference_mask}")
reference_json = f"{reference_prefix}.json"

#%%


brite_urls = [
    i for i in parse_links_from_soup(get_soup("https://www.genome.jp/kegg/brite.html"))
    if re.match("/brite/ko[0-9]+$", i) is not None
]
ko_ids = sorted(set([re.findall(KO_REGEX, i)[0] for i in brite_urls]))
ko_ids

#%%


reference_dicts_list = multi_core_queue2(_get_dict, ko_ids)

#%%

out_dicts = list()
# for reference_dict in reference_dicts_list:
#     _expand_kegg_dict(reference_dict)
reference_dict = [i for i in reference_dicts_list if i.get("name") == "ko00001"][0]
_expand_kegg_dict(reference_dict)

raw_reference_df = pd.DataFrame(out_dicts)
reference_df = raw_reference_df.loc[
    raw_reference_df["ko"].str.startswith("K"),
    :
].set_index("ko").sort_index().drop_duplicates()
reference_df

#%%


normalized_reference_df = reference_df.fillna("").groupby(reference_df.index).agg(
    lambda x: " && ".join(sorted(set(remove_empty_values(x))))
)
normalized_reference_df

#%%


qiime_styled_column_mapper = list(TAXONOMY_CHARS.keys())[:normalized_reference_df.shape[1]]
qiime_styled_df_raw = normalized_reference_df.rename_axis(index="#OTU ID").rename(
    columns={
        k: v
        for k, v
        in zip(normalized_reference_df.columns, qiime_styled_column_mapper)
    }
)
qiime_styled_df_raw

#%%


qiime_styled_df = pd.DataFrame(
    qiime_styled_df_raw.apply(
        lambda row: '; '.join([
            f"{column_name}__{value}"
            for column_name, value in row.iteritems()
        ]),
        axis=1
    ).rename("taxonomy")
)
qiime_styled_df

#%%


dump_dict(reference_dicts_list, reference_json)
dump_tsv(reference_df, f"{reference_prefix}_denormalized.tsv", reset_index=True)
dump_tsv(normalized_reference_df, f"{reference_prefix}_normalized.tsv", reset_index=True)
dump_tsv(qiime_styled_df, f"{reference_prefix}_qiime_styled.tsv", reset_index=True)
