import polars as pl

from typing import TextIO
from censtats.length.constants import (
    DEF_INPUT_BED_COLS,
    DEF_INPUT_RM_COLS,
    DEF_INPUT_RM_COL_IDX,
)


def read_rm(infile: TextIO | str) -> pl.DataFrame:
    """
    Read a tab-delimited RepeatMasker output file and adjust relative coordinates to absolute coordinates.
    """
    return (
        pl.read_csv(
            infile,
            separator="\t",
            has_header=False,
            columns=DEF_INPUT_RM_COL_IDX,
            new_columns=DEF_INPUT_RM_COLS,
            truncate_ragged_lines=True,
        )
        .with_columns(
            ctg_name=pl.col("contig").str.extract(r"^(.*?):|^(.*?)$"),
            ctg_st=pl.col("contig").str.extract(r":(\d+)-").cast(pl.Int64).fill_null(0),
            ctg_end=pl.col("contig")
            .str.extract(r"-(\d+)$")
            .cast(pl.Int64)
            .fill_null(0),
        )
        # Adjust for contig coordinates if any.
        .with_columns(
            pl.col("start") + pl.col("ctg_st"), pl.col("end") + pl.col("ctg_st")
        )
    )


def read_stv(infile: TextIO | str) -> pl.DataFrame:
    """
    Read an HOR Stv bed.
    """
    return pl.read_csv(
        infile,
        separator="\t",
        columns=[0, 1, 2, 3, 4, 5],
        new_columns=DEF_INPUT_BED_COLS,
        has_header=False,
    )


def format_and_output_lengths(
    df: pl.DataFrame,
    output: TextIO | str,
    output_cols: list[str],
) -> None:
    (
        df.with_columns(
            sort_idx=pl.col("chrom")
            .str.extract(r"chr([0-9XY]+)")
            .replace({"X": "23", "Y": "24"})
            .cast(pl.Int32)
        )
        .sort(by="sort_idx")
        .select(output_cols)
        .write_csv(output, include_header=False, separator="\t")
    )
