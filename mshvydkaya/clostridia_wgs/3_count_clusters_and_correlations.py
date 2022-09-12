#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%

import os
import pandas as pd
from matplotlib import pyplot as plt
from meta.utils.io import dump_string
from meta.utils.diversity import get_newick
from meta.utils.primitive import flatten_2d_array
from meta.utils.pandas import dump_tsv, excel_to_dfs_dict
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from meta.utils.correlations import slice_and_correlate_df

#%%

INDEX = "sample_name"
PROJECT_DIR = "/data03/bio/projects/mshvydkaya/clostridia_wgs_13"

coverage_dfs = excel_to_dfs_dict(os.path.join(PROJECT_DIR, "pga-pe-pipeline/coverages.xlsx"))
phenotype_dfs = excel_to_dfs_dict(os.path.join(PROJECT_DIR, "gabrichevsky tables.xlsx"))

#%%

def retrieve_some_db(s: str):
    return coverage_dfs[[i for i in coverage_dfs.keys() if i.startswith(s)][0]]

#%%

tables_by_sample_name = dict()

tables_by_sample_name["card"] = retrieve_some_db("card").loc[
    :,
    [
        INDEX,
        "Best_Hit_ARO",
        "Best_Identities"
    ]
].pivot_table(
    index=INDEX,
    columns="Best_Hit_ARO",
    values="Best_Identities"
)

tables_by_sample_name["mvirdb"] = retrieve_some_db("mvirdb").loc[
    :,
    [
        INDEX,
        "feature_names",
        "id_mapped_relative_abundance"
    ]
].pivot_table(
    index=INDEX,
    columns="feature_names",
    values="id_mapped_relative_abundance"
)

tables_by_sample_name["tadb"] = retrieve_some_db("tadb").loc[
    :,
    [
        INDEX,
        "protein_description",
        "id_mapped_relative_abundance"
    ]
].pivot_table(
    index=INDEX,
    columns="protein_description",
    values="id_mapped_relative_abundance"
)

tables_by_sample_name["vfdb"] = retrieve_some_db("vfdb").loc[
    :,
    [
        INDEX,
        "protein_description",
        "id_mapped_relative_abundance"
    ]
].pivot_table(
    index=INDEX,
    columns="protein_description",
    values="id_mapped_relative_abundance"
)

#%%

sample_names = sorted(
    set(
        flatten_2d_array([
            i.index.values for i in tables_by_sample_name.values()
        ])
    )
)

#%%

phenotype_dfs = {
    k: v for k, v in zip(
        ["common_results", "drug_test_results"],
        phenotype_dfs.values()
    )
}

tables_by_sample_name["common_results"] = pd.get_dummies(
    phenotype_dfs["common_results"].set_index(INDEX)
)

#%%

tables_by_sample_name["drug_test_results_zone_mm"] = phenotype_dfs["drug_test_results"].loc[
    :,
    [
        INDEX,
        "antibiotic",
        "zone_mm"
    ]
].pivot_table(
    index=INDEX,
    columns="antibiotic",
    values="zone_mm"
)

# Not including `reference_mm`

phenotype_dfs["drug_test_results"]["is_resisted"] = (
        phenotype_dfs["drug_test_results"]["interpretation"] == "resistance"
).astype(int)

tables_by_sample_name["drug_test_results_is_resisted"] = phenotype_dfs["drug_test_results"].loc[
    :,
    [
        INDEX,
        "antibiotic",
        "is_resisted"
    ]
].pivot_table(
    index=INDEX,
    columns="antibiotic",
    values="is_resisted"
)

#%%

tables_by_sample_name = {
    k: v.rename(columns={i: f"{k}@{i}" for i in v.columns})
    for k, v in tables_by_sample_name.items()
}
main_pivot_df = pd.concat(tables_by_sample_name.values(), axis=1).fillna(0).loc[sample_names, :]

#%%

dendrogram_dir = os.path.join(PROJECT_DIR, "dendrogram")
dump_tsv(
    main_pivot_df,
    os.path.join(dendrogram_dir, "main_pivot_table.tsv"),
    reset_index=True,
)

#%%

linkage_data = linkage(
    main_pivot_df,
    method="ward",
    metric="euclidean"
)

tree = to_tree(linkage_data, False)
newick_string = get_newick(tree, tree.dist, leaf_names=main_pivot_df.index)
dump_string(newick_string, os.path.join(dendrogram_dir, "dendrogram.newick"))


#%%

plt.clf()
plt.close()
plt.rcParams.update({
    "figure.figsize": (10, 10),
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})
fig, ax = plt.subplots()


title = f"Dendrogram of experimental data"

_ = dendrogram(
    linkage_data,
    ax=ax,
    labels=main_pivot_df.index,
    orientation="right"
)
ax.set_xlabel("Distance")
ax.set_ylabel(INDEX)
ax.set_title(title)

plt.savefig(os.path.join(dendrogram_dir, "dendrogram.pdf"))
plt.clf()
plt.close()

#%%

correlation_dfs = slice_and_correlate_df(main_pivot_df)

#%%

correlation_dir = os.path.join(PROJECT_DIR, "correlation")

_ = [
    dump_tsv(v, os.path.join(correlation_dir, f"{k}.tsv"))
    for k, v in correlation_dfs.items()
]
