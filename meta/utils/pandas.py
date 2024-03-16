#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import traceback
import numpy as np
import pandas as pd
import joblib as jb
from scipy import stats
from typing import (
    Callable,
    Iterable,
    List,
)


def load_tsv(table, **kwargs):
    _kwargs = dict(encoding="utf-8", engine="python", header=0, sep="\t")
    if len(kwargs.keys()) > 0:
        _kwargs.update(kwargs)
    return pd.read_csv(table, **_kwargs)


def dump_tsv(
    df: pd.DataFrame,
    table_file: str,
    col_names: list = None,
    reset_index: bool = False
):
    assert isinstance(df, pd.DataFrame)
    _df = df.copy()
    os.makedirs(os.path.dirname(table_file), exist_ok=True)
    if col_names is not None and len(col_names) > 0:
        _df = _df.loc[:, col_names]
    if reset_index:
        _df.reset_index(inplace=True)
    _df.to_csv(table_file, encoding="utf-8", sep="\t", index=False, header=True)


def dict2pd_series(
    dictionary: dict,
    dtype: str = "object",
    name: str = "",
    sort_keys: bool = False
):
    out = pd.Series(dtype=dtype)
    keys = list(dictionary.keys())
    if sort_keys:
        keys = sorted(keys)
    for key in keys:
        out.at[key] = dictionary[key]
    if len(name) > 0:
        out.name = name
    return out


def apply_mp_function_to_df(
    func: Callable,
    df: pd.DataFrame,
    index_name: str = "",
    columns_name: str = "",
    **kwargs
):
    if len(index_name) == 0:
        index_name = df.columns.name
    if len(columns_name) == 0:
        columns_name = df.index.name
    results = jb.Parallel(n_jobs=-1)(
        jb.delayed(func)(j, **kwargs)
        for j in [df[i] for i in df.columns]
    )
    out_df = pd.DataFrame(results)
    out_df = out_df.reindex(
        sorted(out_df.columns), axis=1
    ).rename_axis(
        index=index_name, columns=columns_name
    ).transpose()
    return out_df


def apply_mp_pairwise_function_to_df(
    func: Callable,
    df: pd.DataFrame,
    columns_name: str = "",
    value_column_name: str = "value",
    **kwargs
):
    from itertools import combinations
    from warnings import catch_warnings, filterwarnings
    from meta.utils.diversity import wrapper_for_pairwise_function
    from meta.utils.queue import multi_core_queue2

    def _mp_apply_pairwise_function(combination: List[str]):
        # Distance functions may produce a lot of warnings
        with catch_warnings():
            filterwarnings("ignore")
            value = wrapper_for_pairwise_function(
                func=func,
                x=df.loc[combination, :].values.tolist(),
                **kwargs,
            )
        d = {
            "x": combination[0],
            "y": combination[-1],
            value_column_name: value
        }
        return dict2pd_series(d, dtype="float64")

    if len(columns_name) == 0:
        columns_name = df.index.name
    pairwise_out = multi_core_queue2(
        _mp_apply_pairwise_function,
        combinations(df.index, 2)
    )
    out_df = pd.DataFrame(pairwise_out).set_index(
        ["x", "y"]
    ).rename_axis(columns=columns_name)
    return out_df


def normalize_series(s: pd.Series):
    return 100 * s / s.sum()


def normalize_df(df: pd.DataFrame):
    return apply_mp_function_to_df(func=normalize_series, df=df)


def corr(df: pd.DataFrame, methods: list = None):
    if methods is None or len(methods) != 2:
        methods = [lambda x, y: stats.spearmanr(x, y)[0], lambda x, y: stats.spearmanr(x, y)[-1]]
    correlation_df = df.corr(method=methods[0])
    p_values_df = df.corr(method=methods[-1]) - np.eye(*correlation_df.shape)
    asterisk_df = p_values_df.applymap(
        lambda x: "".join(["*" for i in (0.01, 0.05, 0.1) if x < i]))
    denoted_correlation_df = correlation_df.round(2).astype(str) + asterisk_df
    d = {"correlations": correlation_df, "p_values": p_values_df,
         "denoted_correlation": denoted_correlation_df}
    try:
        d["name"] = df.name
    except AttributeError:
        pass
    return d


def convert_objects(df: pd.DataFrame, type_=str):
    _df = df.copy()
    for column in _df.columns:
        _df[column] = _df[column].astype(type_)
    return _df


def find_duplicated_rows(df: pd.DataFrame, on: str):
    return df.loc[df[on].duplicated(keep=False)].sort_values(on)


def find_notna_rows(df: pd.DataFrame, on: str):
    return df.loc[~df[on].isna()].sort_values(on)


def deduplicate_df_by_row_merging(df: pd.DataFrame, on: str, sep: str = ";"):
    def _deduplicate(series: pd.Series):
        return sep.join(sorted(set([str(i) for i in series.values])))
    _df = convert_objects(df)
    duplicated_df = find_duplicated_rows(_df, on)
    out_series = []
    for duplicated_value in sorted(set(duplicated_df[on].values)):
        out_series.append(
            duplicated_df.loc[duplicated_df[on] == duplicated_value].apply(_deduplicate, axis=0)
        )
    return pd.concat(
        [pd.DataFrame(out_series), _df.drop(index=duplicated_df.index)],
        axis=0
    ).sort_values(on)


def concat(
    dfs: Iterable[pd.DataFrame],
    on,
    axis: int = 0,
    how: str = "outer",
    columns_name: str = "",
    set_index: bool = True,
    reset_index: bool = True,
    fillna=None,
    **kwargs
):
    assert all(isinstance(i, pd.DataFrame) for i in dfs)
    _dfs = list(dfs)
    if set_index:
        _dfs = [i.set_index(on) for i in _dfs]
    out = pd.concat(
        _dfs, join=how, axis=axis, **kwargs
    ).rename_axis(index=on, columns=columns_name).sort_index()
    if fillna is not None:
        out = out.fillna(fillna)
    if reset_index:
        out = out.reset_index()
    return out


def merge(left_df: pd.DataFrame, right_df: pd.DataFrame, on: str, how: str = "outer",
          update: bool = True, deduplicate: bool = False, **kwarg):
    assert all(isinstance(i, pd.DataFrame) for i in [left_df, right_df])
    _left_df = left_df.copy()
    _right_df = right_df.copy()
    if update:
        _left_df.update(_right_df)
    df1_unique_columns = [i for i in _right_df.columns if i not in _left_df.columns]
    out = _left_df.merge(_right_df.loc[:, [on] + df1_unique_columns], on=on, how=how, **kwarg)
    if deduplicate:
        out = deduplicate_df_by_row_merging(out, on)
    return out


def count_column_sizes(df: pd.DataFrame):
    return df.apply(lambda y: max(y.map(lambda x: len(str(x)))))


def excel_to_dfs_dict(file: str, **kwargs):
    sheets = pd.ExcelFile(file).sheet_names
    out = dict()
    for sheet in sheets:
        try:
            out[sheet] = pd.read_excel(file, sheet_name=sheet, **kwargs)
        except Exception:
            traceback.print_exc()
            print(f"Cannot parse sheet '{sheet}' from file '{file}'")
    return out


def dfs_dict_to_excel(d: dict, file: str, **kwargs):
    """
    :param d: {sheet_name <str>: sheet_dataframe <pandas.DataFrame>}
    :param file: table.xlx
    :return:
    """
    import xlsxwriter
    w = pd.ExcelWriter(file, engine="xlsxwriter")
    for sheet, df in d.items():
        try:
            df.to_excel(w, index=False, sheet_name=sheet, **kwargs)
        except Exception:
            traceback.print_exc()
            print(f"Cannot save sheet '{sheet}' to file '{file}'")
    w.save()


def remove_longest_columns(df: pd.DataFrame, size: int = 32767):  # M$ Excel cell size limit
    max_lengths = count_column_sizes(df)
    return df.loc[:, max_lengths.loc[max_lengths < size].index]


def dwell_df_on_column(df: pd.DataFrame, column_name: str):
    if column_name not in df.columns:
        return df
    return df.loc[:, [column_name] + [i for i in df.columns if i != column_name]]


def get_major_features_df(
    df: pd.DataFrame,
    n: int = 10,
    other_column_name: str = "Others"
):
    # Input df's columns: features, indexes: samples
    sum_series = df.sum().sort_values(ascending=False)
    majors = sum_series[:n].index
    others = sum_series[n:].index
    out_df = df.loc[:, majors]
    out_df[other_column_name] = df[others].sum(axis=1)
    return out_df


def split_df_into_chunks_of_size(
    df: pd.DataFrame,
    axis: int = 0,
    chunk_size: int = 10,
    separator: str = "_to_",
):
    from meta.utils.primitive import split_list_into_chunks_of_size
    if axis >= 2:
        raise ValueError(f"Axis must be 0 or 1 ({axis} is given)")
    index_lists = split_list_into_chunks_of_size(
        list(range(df.shape[axis])), chunk_size
    )
    out = dict()
    for index_list in index_lists:
        first = index_list[0]
        last = index_list[-1]
        # .iloc[] is primarily integer position based (from 0 to length-1 of the axis)
        if axis == 0:
            df_1 = df.iloc[first:last + 1, :]
        else:
            df_1 = df.iloc[:, first:last + 1]
        if df_1.shape[axis] == 0:
            continue
        df_1.reindex(
            sorted(df_1.axes[axis]), axis=axis
        )
        name = f"{df_1.axes[axis][0]}{separator}{df_1.axes[axis][-1]}"
        df_1.name = name
        out[name] = df_1
    return out


def draw_supervenn_diagram(
    data: dict,
    title: str,
    output_dir: str,
    diagram_labels: Iterable = None
):
    """
    :param data: dict,
        keys are categories (e.g. sample names, sources, etc.)
        values are sets of non-empty feature names
    :param title:
    :param output_dir:
    :param x_label:
    :param y_label:
    :return:
    """
    import seaborn as sns
    from supervenn import supervenn
    from matplotlib import pyplot as plt
    from meta.utils.io import dump_dict
    plt.clf()
    plt.close()

    sns.set(style="whitegrid")
    sns.set_palette(os.getenv("MATPLOTLIB_COLORMAP", "hsv"))

    plt.rcParams.update({
        "figure.figsize": (5, 5),
        "figure.dpi": 75
    })

    supervenn(
        sets=data.values(),
        set_annotations=list(data.keys()),
        side_plots=False,
        widths_minmax_ratio=0.05,
    )
    if type(diagram_labels) in (tuple, list) and len(diagram_labels) > 1:
        plt.xlabel(diagram_labels[0])
        plt.ylabel(diagram_labels[1])

    plt.suptitle(title)
    # plt.show()
    file_mask = os.path.join(output_dir, title)
    os.makedirs(output_dir, exist_ok=True, mode=0o777)
    plt.tight_layout()
    plt.savefig(f"{file_mask}_supervenn_diagram.jpg")
    plt.clf()
    plt.close()
    dump_dict({k: list(v) for k, v in data.items()}, f"{file_mask}_data.json")
    return True


def count_features_and_draw_supervenn_diagram(
    df: pd.DataFrame,
    title: str,
    output_dir: str,
    diagram_labels: None
):
    """
    :param df:
        indexes are categories (e.g. sample names, sources, etc.)
        columns are features
    :param title:
    :param output_dir:
    :param diagram_labels:
    :return:
    """
    from meta.utils.primitive import remove_empty_values

    df = df.loc[:, df.sum() > 0]

    row_count_df = df.apply(
        lambda x: x.apply(
            lambda y: x.name if y > 0 else ""
        )
    ).astype(str)

    row_count_dict = {
        k: set(remove_empty_values(v.values))
        for k, v in row_count_df.transpose().to_dict("series").items()
    }

    column_count_dicts = list()
    for column_name in row_count_df:
        feature_occurrences = [
            str(k) for k, v in row_count_dict.items() if column_name in v
        ]
        column_count_dicts.append(dict(
            feature=column_name,
            occurrences=",".join(feature_occurrences),
            frequency=len(feature_occurrences),
        ))
    column_count_df = pd.DataFrame(column_count_dicts).sort_values(
        "frequency", ascending=False
    ).set_index("feature")
    column_count_df = pd.concat([
        column_count_df,
        df.transpose(),
    ], axis=1, sort=False)
    column_count_df.rename_axis(
        index=df.columns.name,
        columns=df.index.name,
        inplace=True
    )
    column_count_df.name = title

    file_mask = os.path.join(output_dir, title)
    dump_tsv(
        column_count_df,
        f"{file_mask}_frequencies.tsv",
        reset_index=True
    )

    draw_supervenn_diagram(
        data=row_count_dict,
        title=title,
        output_dir=output_dir,
        diagram_labels=diagram_labels,
    )
    return column_count_df


def grouping_describe(df: pd.DataFrame):
    """
    :param df:
        indexes are categories (non-unique are preferred by design)
        columns are features
    :return:
    """
    def _process(s: str):
        return df.loc[s, :].describe().rename(
            axis=0,
            mapper=lambda x: f"{s}-{x}"
        ).transpose()

    dfs = jb.Parallel(n_jobs=-1)(
        jb.delayed(_process)(i)
        for i in sorted(set(df.index))
    )

    out = pd.concat(dfs, axis=1, sort=False)
    return out


def count_feature_based_group_relations(
    df: pd.DataFrame,
    grouping_column_name: str,
    feature_name: str,
    output_dir: str,
    annotation_df: pd.DataFrame = None
):
    """
    :param df:
        indexes are categories (e.g. sample names, sources, etc.)
        columns are features AND 1 grouping column
    :param grouping_column_name:
    :param feature_name:
    :param output_dir:
    :param annotation_df:
    :return:
    """
    def _process(rounding_policy: str, feature_grouped_df: pd.DataFrame):
        out_dict = dict()
        value_count_dir = os.path.join(
            output_dir,
            f"{feature_name}_value_counts"
        )
        unannotated_feature_frequency_df = count_features_and_draw_supervenn_diagram(
            df=feature_grouped_df,
            title=f"{feature_name} occurring {rounding_policy} for each sample per sample group",
            output_dir=value_count_dir,
            diagram_labels=(f"Distinct {feature_name}", grouping_column_name),
        )
        out_dict["unannotated_feature_frequency_df"] = unannotated_feature_frequency_df
        if isinstance(annotation_df, pd.DataFrame) and annotation_df.shape[0] > 0:
            annotated_feature_frequency_df = pd.concat(
                [
                    unannotated_feature_frequency_df,
                    annotation_df
                ],
                axis=1,
                join="inner",
                sort=False,
            ).sort_values("frequency", ascending=False)
            dump_tsv(
                annotated_feature_frequency_df,
                os.path.join(
                    value_count_dir,
                    f"{unannotated_feature_frequency_df.name}_frequencies_annotated.tsv"
                ),
                reset_index=True
            )
            out_dict["annotated_feature_frequency_df"]: annotated_feature_frequency_df
        return out_dict

    if not isinstance(annotation_df, pd.DataFrame):
        annotation_df = pd.DataFrame()

    annotation_df = pd.concat(
        [
            grouping_describe(df.set_index(grouping_column_name, append=False)),
            annotation_df,
        ],
        axis=1,
        join="inner",
        sort=False,
        names=[feature_name]
    )

    feature_value_count_series = df[grouping_column_name].value_counts()

    feature_grouped_df_dict = dict()
    feature_counter_df = (
        df.set_index(grouping_column_name, append=True) > 0
    ).astype(int).reset_index(grouping_column_name)
    count_feature_grouped_df = feature_counter_df.groupby(grouping_column_name).sum()
    percentage_feature_grouped_df = count_feature_grouped_df.divide(
        feature_value_count_series, axis=0
    ).rename_axis(index=grouping_column_name, columns=df.index.name)

    feature_grouped_df_dict["at least once"] = percentage_feature_grouped_df
    feature_grouped_df_dict["everywhere"] = percentage_feature_grouped_df.astype(int)

    d = jb.Parallel(n_jobs=-1)(
        jb.delayed(_process)(k, v)
        for k, v in feature_grouped_df_dict.items()
    )
    return d


def draw_pca(
    df: pd.DataFrame,
    name: str,
    grouping_column_name: str,
    output_dir: str,
):
    """
    :param df:
        indexes are categories (e.g. sample names, sources, etc.)
        columns are features AND 1 grouping column
    :param name:
    :param grouping_column_name:
    :param output_dir:
    :return:
    """
    from matplotlib import pyplot as plt
    from meta.utils.io import dump_dict
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    import seaborn as sns

    pca_components = 2
    scaling_df = df.set_index(grouping_column_name, append=True)

    scaler = StandardScaler()
    scaled_array = scaler.fit_transform(scaling_df)

    pca = PCA(n_components=pca_components)
    pca_array = pca.fit_transform(scaled_array)

    pca_column_names = [f"PCA {i + 1}" for i in range(pca_array.shape[1])]

    pca_df = pd.DataFrame(
        pca_array,
        index=scaling_df.index,
        columns=pca_column_names
    ).reset_index(grouping_column_name)

    plt.clf()
    plt.close()

    sns.set(style="whitegrid")

    sns.set_palette(os.getenv("MATPLOTLIB_COLORMAP", "turbo"))

    plt.rcParams.update({
        "figure.figsize": (5, 5),
        "figure.dpi": 75
    })

    sns.lmplot(
        **{k: v for k, v in zip(["x", "y"], pca_column_names)},
        data=pca_df,
        hue=grouping_column_name,
        fit_reg=False,
        legend=True
    )

    nd = len(pca_column_names)
    cd = len(set(df[grouping_column_name]))

    title = f"{nd}D PCA for {name}  ({cd} groups)"
    plt.title(title)

    file_mask = os.path.join(output_dir, title)
    dump_tsv(pca_df, f"{file_mask}.tsv", reset_index=True)
    dump_dict(scaled_array.tolist(), f"{file_mask}_scaled_values.json")
    plt.savefig(f"{file_mask}.jpg")

    plt.clf()
    plt.close()


def draw_pcoa(
    df: pd.DataFrame,
    name: str,
    grouping_column_name: str,
    output_dir: str,
):
    """
    :param df:
        indexes are categories (e.g. sample names, sources, etc.)
        columns are features AND 1 grouping column
    :param name:
    :param grouping_column_name:
    :param output_dir:
    :return:
    """
    from matplotlib import pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import normalize
    import seaborn as sns
    from meta.utils.io import dump_dict

    pca_components = 2
    scaling_df = df.set_index(grouping_column_name, append=True)

    normalized_array = normalize(scaling_df).transpose()

    pca = PCA(n_components=pca_components)
    pca_array = pca.fit(normalized_array)
    pcoa_array = pca.components_.transpose()

    pcoa_column_names = [f"PCoA {i + 1}" for i in range(pcoa_array.shape[1])]

    pcoa_df = pd.DataFrame(
        pcoa_array,
        index=scaling_df.index,
        columns=pcoa_column_names
    ).reset_index(grouping_column_name)

    plt.clf()
    plt.close()

    sns.set(style="whitegrid")

    sns.set_palette(os.getenv("MATPLOTLIB_COLORMAP", "turbo"))

    plt.rcParams.update({
        "figure.figsize": (5, 5),
        "figure.dpi": 75
    })

    sns.lmplot(
        **{k: v for k, v in zip(["x", "y"], pcoa_column_names)},
        data=pcoa_df,
        hue=grouping_column_name,
        fit_reg=False,
        legend=True
    )

    nd = len(pcoa_column_names)
    cd = len(set(df[grouping_column_name]))

    title = f"{nd}D PCoA for {name}  ({cd} groups)"
    plt.title(title)

    file_mask = os.path.join(output_dir, title)
    dump_tsv(pcoa_df, f"{file_mask}.tsv", reset_index=True)
    dump_dict(normalized_array.tolist(), f"{file_mask}_normalized_values.json")
    plt.savefig(f"{file_mask}.jpg")

    plt.clf()
    plt.close()


PRINCIPAL_ANALYSIS_METHODS = dict(draw_pca=draw_pca, draw_pcoa=draw_pcoa,)


def annotate_and_aggregate_df(df: pd.DataFrame, annotation_df: pd.DataFrame):
    """
    :param df:
        indexes are features
        columns are categories (e.g. sample names, sources, etc.)
    :param annotation_df:
        indexes are features (non-unique are preferred by design)
        columns are grouping / aggregation columns
    :return:
    """
    index_column_name = df.index.name
    annotated_df = annotation_df.reset_index().merge(
        df.reset_index(), how="right", on=index_column_name
    ).set_index(index_column_name)
    value_column_names = df.columns.tolist()
    out_dict = dict()
    for grouping_column_name in annotation_df.columns:
        grouped_df = annotated_df.loc[
            :,
            [grouping_column_name] + value_column_names
        ].groupby(grouping_column_name).sum()
        grouped_df = pd.concat(
            [
                annotated_df.loc[
                    :,
                    [grouping_column_name]
                ].reset_index().groupby(grouping_column_name).apply(
                    lambda x: ";".join(str(i) for i in sorted(set(x[index_column_name].values)))
                ).rename(index_column_name),
                grouped_df,
            ], axis=1, sort=False
        )
        out_dict[grouping_column_name] = grouped_df
    return annotated_df, out_dict
