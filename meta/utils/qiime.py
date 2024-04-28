
import pandas as pd


def fix_sample_ids(df: pd.DataFrame):
    # Deblur cannot operate on sample IDs that contain underscores
    df["#SampleID"].replace(
        "[^A-Za-z0-9]+",
        "-",
        inplace=True,
        regex=True,
    )
