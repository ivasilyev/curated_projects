#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import joblib as jb
import pandas as pd
from numpy import nan
from skbio import TreeNode
from skbio.diversity import alpha as a
from skbio.diversity import beta as b
from scipy.spatial import distance as d
from meta.utils.pandas import apply_mp_function_to_df, dict2pd_series


TAXONOMY_ORDER = "Domain, Phylum, Class, Order, Family, Genus, Species".split(", ")
taxonomy_chars = {i[0].lower(): i for i in TAXONOMY_ORDER}
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


def parse_taxonomy(s: str):
    # s = "d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia-Shigella; s__Escherichia_coli"
    out = {j[0].strip(): "__".join(j[1:]).strip() for j in [i.split("__") for i in s.split("; ")]}
    out = {k if k not in taxonomy_chars.keys() else taxonomy_chars[k]: out[k] for k in out}
    return out


def root_tree_node(tree: TreeNode):
    import warnings
    from skbio.util import RepresentationWarning
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RepresentationWarning)
        return tree.root_at_midpoint()


def count_alpha_diversity(series: pd.Series, rooted_tree: TreeNode = None):
    d = {k: ALPHA_DIVERSITY_FUNCTIONS[k](series.values) for k in ALPHA_DIVERSITY_FUNCTIONS}
    d["Inverse Simpson Index"] = 1.0 / d["Simpson Index"]
    d["Gini–Simpson Index"] = 1.0 - d["Simpson Index"]
    if rooted_tree is not None:
        d["Faith Diversity"] = a.faith_pd(series.values, series.index.values, rooted_tree)
    out = dict2pd_series(d)
    out.name = series.name
    return out


def count_alpha_diversity_df(df: pd.DataFrame, rooted_tree: TreeNode = None):
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
