import polars as pl
from typing import TextIO
from .constants import RM_COLS, RM_COL_IDX


def read_repeatmasker_output(input_path: TextIO | str) -> pl.DataFrame:
    return pl.read_csv(
        input_path,
        separator="\t",
        has_header=False,
        columns=RM_COL_IDX,
        new_columns=RM_COLS,
        truncate_ragged_lines=True,
    ).with_columns(dst=pl.col("end") - pl.col("start"))
