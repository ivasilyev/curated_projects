
import os
import joblib as jb
import pandas as pd
from meta.utils.pandas import dump_tsv, load_tsv
from meta.utils.diversity import parse_taxonomy
from meta.utils.file_system import find_file_by_tail

reference_db_name = "SILVA"
output_dir = "/data/reference"

reference_version = 138
reference_mask = f"{reference_db_name}_v{reference_version}"
reference_dir = os.path.join(output_dir, reference_db_name, reference_mask)
reference_prefix = os.path.join(reference_dir, reference_mask)

features_df = load_tsv(find_file_by_tail(reference_dir, "Taxonomy_headed.tsv"))


def _expand(d: dict):
    # Keys: #OTU ID, taxonomy
    return dict(**i for i in [d, parse_taxonomy(d["taxonomy"])])


taxa_dicts = jb.Parallel(n_jobs=-1)(
    jb.delayed(_expand)(i)
    for i in features_df.to_dict("records")
)


expanded_df = pd.DataFrame(taxa_dicts).sort_values("#OTU ID")
dump_tsv(expanded_df, f"{reference_prefix}_taxonomy_expanded.tsv")
