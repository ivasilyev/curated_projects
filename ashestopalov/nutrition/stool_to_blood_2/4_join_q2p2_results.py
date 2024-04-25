import os
from typing import Dict

import pandas as pd

from meta.utils.file_system import find_file_by_tail
from meta.utils.pandas import df_to_7z, concat, load_tsv
from meta.utils.primitive import dicts_list_to_lists_dict
from meta.utils.queue import multi_core_queue2


def join_q2p2_results(split_results_dir: str):
    def _mp_join_q2p2_result(table_file: str) -> Dict[str, pd.DataFrame]:
        sample_column_name = "#SampleID"
        if table_file == "EC_metagenome_out_pred_metagenome_unstrat_described.tsv":
            index_columns = ["function", "description"]
            df = concat(
                dfs=split_results_tables_by_basename.get(table_file),
                on=index_columns,
                axis=1,
                columns_name=sample_column_name,
                fillna=0,
            ).sort_values(index_columns[-1]).set_index(index_columns)
            return {"picrust2_ec": df}
        elif table_file == "KO_metagenome_out_pred_metagenome_unstrat_described.tsv":
            index_columns = ["function", "description"]
            df = concat(
                dfs=split_results_tables_by_basename.get(table_file),
                on=index_columns,
                axis=1,
                columns_name=sample_column_name,
                fillna=0,
            ).sort_values(index_columns[-1]).set_index(index_columns)
            return {"picrust2_ko": df}
        elif table_file == "path_abun_unstrat_described.tsv":
            index_columns = ["pathway", "description"]
            df = concat(
                dfs=split_results_tables_by_basename.get(table_file),
                on=index_columns,
                axis=1,
                columns_name=sample_column_name,
                fillna=0,
            ).sort_values(index_columns[-1]).set_index(index_columns)
            return {"picrust2_pathway": df}
        elif table_file == "OTUs_with_taxa.tsv":
            index_columns = ["#OTU ID", "taxonomy"]
            df = concat(
                dfs=split_results_tables_by_basename.get(table_file),
                on=index_columns,
                axis=1,
                fillna=0,
            ).sort_values(index_columns[-1]).set_index(index_columns)
            return {"qiime2_otu": df}

    # %%

    split_results_tables = [
        i for i in
        find_file_by_tail(split_results_dir, ".tsv", True)
        if os.path.basename(os.path.dirname(i)) == "results"
           and not i.endswith("pred_metagenome_contrib.legacy.tsv")
    ]

    # split_results_tables

    # %%

    split_results_dfs = multi_core_queue2(
        lambda x: {os.path.basename(x): load_tsv(x).sort_index(axis=1)},
        split_results_tables
    )

    # split_results_dfs

    # %%

    split_results_tables_by_basename = dicts_list_to_lists_dict(split_results_dfs)
    split_results_tables_by_basename.keys()

    joined_dfs = multi_core_queue2(
        _mp_join_q2p2_result,
        split_results_tables_by_basename.keys(),
    )

    joined_dfs_dict = dicts_list_to_lists_dict(joined_dfs, truncate=0)
    joined_dfs_dict.keys()

    pipeline_out_dir = os.path.join(
        split_results_dir,
        "qiime2-picrust2-pipeline",
        "pipeline_out"
    )
    queue = [
        dict(
            df=v.reset_index(),
            archive=os.path.join(pipeline_out_dir, f"{k}.7z"),
            table_basename=k,
        )
        for k, v in joined_dfs_dict.items()
    ]

    # %%

    # queue
    # queue[0]["df"].columns.tolist()

    # %%

    _ = multi_core_queue2(lambda x: df_to_7z(**x), queue)


join_q2p2_results("/data03/bio/projects/ashestopalov/nutrition/stool_to_blood_2")
