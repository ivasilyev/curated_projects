#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import joblib as jb
import pandas as pd
from numpy import nan
from skbio import TreeNode
from skbio.diversity import alpha as a
from skbio.diversity import beta as b
from scipy.spatial import distance as d


TAXONOMY_RANKS = ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")


def taxa_rank_to_char(s: str):
    s = s.strip(" _")
    if s not in TAXONOMY_RANKS:
        raise ValueError(f"Not a taxonomy_rank: '{s}'")
    return s[0].lower()


TAXONOMY_CHARS = {taxa_rank_to_char(i): i for i in TAXONOMY_RANKS}

ALPHA_DIVERSITY_FUNCTIONS = {
    "Distinct OTUs": a.observed_otus,
    "Shannon Entropy": a.shannon,
    "Berger-Parker Dominance": a.berger_parker_d,
    "Chao1 Richness": a.chao1,
    "Simpson Index": a.simpson
}

BETA_DIVERSITY_FUNCTIONS = {
    "Bray-Curtis Distance": d.braycurtis,
    "Canberra Distance": d.canberra,
    "Chebyshev Distance": d.chebyshev,
    "Manhattan Distance": d.cityblock,
    "Correlation Distance": d.correlation,
    "Cosine Distance": d.cosine,
    "Euclidean Distance": d.euclidean,
    "Jensen-Shannon Metric": d.jensenshannon,
    "Minkowski Distance ": d.minkowski,
    "Squared Euclidean Distance": d.sqeuclidean,
}


def parse_taxonomy(s: str, prefix: str = "__", separator: str = "; " ):
    # s = "d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia-Shigella; s__Escherichia_coli"
    out = {j[0].strip(): prefix.join(j[1:]).strip() for j in [i.split(prefix) for i in s.split(separator)]}
    out = {k if k not in TAXONOMY_CHARS.keys() else TAXONOMY_CHARS[k]: out[k] for k in out}
    return out


def collapse_taxa_df_by_rank(
    df: pd.DataFrame,
    rank: str,
    taxa_column_name: str = "taxonomy",
    is_name_df: bool = True
):
    from meta.utils.pandas import normalize_df
    out_df = normalize_df(
        df.loc[
            :,
            [i for i in df.columns if i != taxa_column_name]
        ].groupby(
            by=df[taxa_column_name].str.extract(
                f"{taxa_rank_to_char(rank)}__([^;]+);"
            )[0].rename(rank)
        ).sum()
    ).dropna(axis=0, how="all")
    if is_name_df:
        out_df.name = rank
    return out_df


def annotate_df_with_taxa_columns(
        taxa_df: pd.DataFrame,
        taxa_column_name: str = "taxonomy",
        columns_to_excude: tuple = ("Unassigned",)
):
    from meta.utils.pandas import series_to_list_of_dicts
    from meta.utils.queue import multi_core_queue2
    def _mp_parse_taxa_columns(d: dict):
        taxonomy_raw = d.get(taxa_column_name)
        taxonomy_parsed_dict = parse_taxonomy(taxonomy_raw)
        return {**d, **taxonomy_parsed_dict}

    taxa_raw_dicts = series_to_list_of_dicts(taxa_df[taxa_column_name])
    taxa_annotated_dicts = multi_core_queue2(_mp_parse_taxa_columns, taxa_raw_dicts)
    annotation_df = pd.DataFrame(taxa_annotated_dicts).set_index(taxa_df.index.name).drop(taxa_column_name, axis=1)
    annotation_df = annotation_df.loc[:, [i for i in annotation_df.columns if i not in columns_to_excude]]
    annotated_taxa_df = pd.concat([taxa_df, annotation_df], axis=1, sort=False)
    return annotated_taxa_df


def root_tree_node(tree: TreeNode):
    import warnings
    from skbio.util import RepresentationWarning
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RepresentationWarning)
        return tree.root_at_midpoint()


def count_alpha_diversity(series: pd.Series, rooted_tree: TreeNode = None):
    from meta.utils.pandas import dict2pd_series
    d = {k: ALPHA_DIVERSITY_FUNCTIONS[k](series.values) for k in ALPHA_DIVERSITY_FUNCTIONS}
    d["Inverse Simpson Index"] = 1.0 / d["Simpson Index"]
    d["Gini–Simpson Index"] = 1.0 - d["Simpson Index"]
    if rooted_tree is not None:
        d["Faith Diversity"] = a.faith_pd(series.values, series.index.values, rooted_tree)
    out = dict2pd_series(d)
    out.name = series.name
    return out


def count_alpha_diversity_df(df: pd.DataFrame, rooted_tree: TreeNode = None):
    """
    :param df:
        indexes are features,
        columns are samples
    :param rooted_tree:
    :return:
    """
    from meta.utils.pandas import apply_mp_function_to_df
    out_df = apply_mp_function_to_df(
        func=count_alpha_diversity,
        df=df,
        index_name=df.columns.name,
        columns_name="alpha_diversity_metric",
        rooted_tree=rooted_tree
    )
    return out_df


def wrapper_for_pairwise_function(func, x: list, **kwargs):
    if len(x) == 2:
        return func(*x, **kwargs)
    else:
        return nan


def count_beta_diversity_df(
    df: pd.DataFrame,  # must contain only sample values & only single `grouping_column_name`
    grouping_column_name: str,
    rooted_tree: TreeNode = None
):
    """
    :param df:
        indexes are samples,
        columns are features
    :param grouping_column_name:
    :param rooted_tree:
    :return:
    """
    def _process(name: str, func, **kwargs):
        _out = dfgb.apply(lambda x: wrapper_for_pairwise_function(
            func=func,
            x=x.drop(grouping_column_name, axis=1).values,
            **kwargs
        )).rename(name)
        return _out
    func_dict = dict(BETA_DIVERSITY_FUNCTIONS)
    dfgb = df.groupby(grouping_column_name)
    metric_dfs = jb.Parallel()(
        jb.delayed(_process)(k, v)
        for k, v in func_dict.items()
    )
    if rooted_tree is not None:
        metric_dfs.append(
            _process("unweighted_unifrac", b.unweighted_unifrac, tree=rooted_tree)
        )
        metric_dfs.append(
            _process("weighted_unifrac", b.weighted_unifrac, tree=rooted_tree)
        )
    out_df = pd.concat(metric_dfs, axis=1, sort=False)
    out_df = out_df.reindex(
        sorted(out_df.columns), axis=1
    ).rename_axis(
        index=grouping_column_name,
        columns=["beta_diversity_metric"]
    )
    return out_df


def get_newick(node, parent_dist, leaf_names, newick="") -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.
    From https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick


def draw_taxa_plots(
    name: str,
    taxa_df: pd.DataFrame,
    taxa_rank: str,
    output_dir: str,
):
    """
    :param name:
    :param taxa_df:
        indexes are samples,
        columns are features
    :param taxa_rank:
    :param output_dir:
    :return:
    """
    import os
    import seaborn as sns
    from matplotlib import pyplot as plt
    from meta.utils.pandas import dump_tsv
    plt.clf()
    plt.close()

    sns.set()

    plt.rcParams.update({
        "figure.figsize": (5, 5),
        "figure.dpi": 75
    })

    subtitle = f"per {taxa_df.shape[0]} samples (from '{name}')"
    ax_1 = taxa_df.plot(
        kind="barh",
        stacked=True,
        title=subtitle,
        mark_right=True,
        colormap="nipy_spectral_r",
        legend=False,
    )

    ax_1.legend(
        loc="upper left",
        bbox_to_anchor=(1.0, 1.0),
        fontsize="small",
    )

    ax_1.set(
        xlabel="Percent in sample",
        ylabel="Sample name",
    )

    title = f"Taxa bar charts collapsed by {taxa_rank}"
    plt.suptitle(title)
    # plt.show()

    output_dir = os.path.join(output_dir, taxa_rank)
    file_mask = os.path.join(
        output_dir, f"{title} {subtitle}."
    )

    os.makedirs(output_dir, exist_ok=True, mode=0o777)
    plt.savefig(f"{file_mask}jpg", bbox_inches="tight")
    dump_tsv(taxa_df, f"{file_mask}tsv", reset_index=True)

    plt.clf()
    plt.close()
    return True


def split_and_draw_taxa_plots(
    collapsed_taxa_df: pd.DataFrame,
    taxa_rank: str,
    chunk_df_size: int,
    output_dir: str,
):
    """
    :param collapsed_taxa_df:
        indexes are samples,
        columns are features
    :param taxa_rank:
    :param chunk_df_size:
    :param output_dir:
    :return:
    """
    from meta.utils.pandas import split_df_into_chunks_of_size
    collapsed_taxa_df_chunks = split_df_into_chunks_of_size(
        collapsed_taxa_df,
        axis=0,
        chunk_size=chunk_df_size,
        separator="' to '"
    )
    for name, taxa_df in collapsed_taxa_df_chunks.items():
        draw_taxa_plots(
            name=name,
            taxa_df=taxa_df,
            taxa_rank=taxa_rank,
            output_dir=output_dir
        )


def collapse_split_and_draw_non_major_df(
    taxa_sample_otu_df: pd.DataFrame,
    major_features_number: int,
    samples_per_time: int,
    output_dir: str
):
    """
    :param taxa_sample_otu_df:
        indexes are features,
        columns are samples AND 1 taxonomy column
    :param major_features_number:
    :param samples_per_time:
    :param output_dir:
    :return:
    """
    from meta.utils.pandas import get_major_features_df

    def _process(taxa_rank: str):
        collapsed_taxa_df = collapse_taxa_df_by_rank(
            df=taxa_sample_otu_df,
            rank=taxa_rank,
            taxa_column_name="taxonomy",
        )

        major_features_df = get_major_features_df(
            df=collapsed_taxa_df.transpose(),
            n=major_features_number,
            other_column_name="Others",
        )

        split_and_draw_taxa_plots(
            collapsed_taxa_df=major_features_df,
            taxa_rank=taxa_rank,
            chunk_df_size=samples_per_time,
            output_dir=output_dir
        )

    _ = jb.Parallel(n_jobs=-1)(
        jb.delayed(_process)(i)
        for i in TAXONOMY_RANKS[:-1]
    )


def draw_dendrogram(
    df: pd.DataFrame,
    name: str,
    output_dir: str,
    metric: str = "braycurtis",
    method: str = "single",
):
    """
    :param df:
        indexes are samples,
        columns are features
    :param name:
    :param output_dir:
    :param metric:
    :param method:
    :return:
    """
    import os
    import seaborn as sns
    from matplotlib import pyplot as plt
    from scipy.spatial.distance import pdist, squareform
    from scipy.cluster.hierarchy import linkage, dendrogram
    from meta.utils.io import dump_dict
    from meta.utils.pandas import dump_tsv

    distance_vector = pdist(
        df.values,
        metric=metric,
    )
    distance_matrix = pd.DataFrame(
        squareform(distance_vector),
        index=df.index,
        columns=df.index,
    )

    tree_df = pd.DataFrame(
        linkage(
            df.values,
            method=method,
            metric=metric
        ),
        columns=["index_1", "index_2", "distance", "new_observations"]
    )

    plt.clf()
    plt.close()
    sns.set()

    plt.rcParams.update({
        "figure.figsize": (20, 5),
        "figure.dpi": 150
    })

    tree_dendrogram = dendrogram(tree_df.values, labels=df.index)

    plt.xlabel("Sample name")
    plt.ylabel("Distance")

    title = f"{name} dendrogram for {df.shape[0]} samples by '{method}'-'{metric}'"
    plt.suptitle(title)
    # plt.show()
    file_mask = os.path.join(output_dir, title)

    os.makedirs(output_dir, exist_ok=True, mode=0o777)
    plt.tight_layout()
    plt.savefig(f"{file_mask}.jpg")
    dump_dict(distance_vector.tolist(), f"{file_mask}_distance_vector.json")
    dump_tsv(distance_matrix, f"{file_mask}_distance_matrix.tsv", reset_index=True)
    dump_tsv(tree_df, f"{file_mask}_tree_table.tsv")
    dump_dict(tree_dendrogram, f"{file_mask}_dendrogram.json")

    plt.clf()
    plt.close()
    return True
