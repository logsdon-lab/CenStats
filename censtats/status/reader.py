import polars as pl

from .constants import RM_COLS


def read_repeatmasker_output(input_path: str) -> pl.LazyFrame:
    return (
        pl.scan_csv(input_path, separator="\t", has_header=False, new_columns=RM_COLS)
        .with_columns(dst=pl.col("end") - pl.col("start"))
        .drop("div", "deldiv", "insdiv", "x", "y", "z", "left", "right", "idx")
    )
