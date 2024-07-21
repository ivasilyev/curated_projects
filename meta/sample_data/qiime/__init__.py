
from pandas import DataFrame
from meta.utils.pandas import load_tsv
from meta.utils.qiime import split_main_metadata_df, dump_sampledata


def split_and_dump_main_sample_data_df(main_metadata_df: DataFrame, sampledata_dir: str):
    split_sampledata_dict = split_main_metadata_df(main_metadata_df)
    dump_sampledata(split_sampledata_dict, sampledata_dir)


def split_and_dump_main_sample_data_table(main_metadata_table: str, sampledata_dir: str):
    main_metadata_df = load_tsv(main_metadata_table)
    split_and_dump_main_sample_data_df(main_metadata_df, sampledata_dir)
