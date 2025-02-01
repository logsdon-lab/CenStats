import sys
import argparse

from typing import TYPE_CHECKING, Any, TextIO

from .estimate_length import hor_array_length
from .constants import (
    DEF_MIN_BLK_HOR_UNITS,
    DEF_MIN_ARR_HOR_UNITS,
    DEF_MIN_ARR_LEN,
    DEF_MIN_ARR_PROP,
    DEF_BP_MERGE_UNITS,
    DEF_BP_MERGE_BLKS,
    DEF_INPUT_BED_COLS,
    DEF_INPUT_RM_COLS,
    DEF_INPUT_RM_COL_IDX,
    DEF_OUTPUT_BED_COLS,
    DEF_OUTPUT_BED_COLS_STRAND,
)
from .io import format_and_output_lengths, read_stv, read_rm

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def add_hor_length_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "length",
        description="Estimate HOR array length from stv bed file / HumAS-HMMER output.",
    )
    ap.add_argument(
        "-i",
        "--input_stv",
        help=f"Input stv row bed file produced by HumAS-HMMER and stv. Expects columns: {DEF_INPUT_BED_COLS}",
        type=argparse.FileType("rb"),
    )
    ap.add_argument(
        "-r",
        "--input_rm",
        help=f"Input tab-delimited RepeatMasker file with no header. Prevents joining across. Expects columns: {DEF_INPUT_RM_COLS} at indices {DEF_INPUT_RM_COL_IDX}.",
        type=argparse.FileType("rb"),
        default=None,
    )
    ap.add_argument(
        "-o",
        "--output",
        help=f"Output bed file with columns: {DEF_OUTPUT_BED_COLS}",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )
    ap.add_argument(
        "-s",
        "--output_strand",
        help=f"Output bed file with columns: {DEF_OUTPUT_BED_COLS}",
        default=None,
        type=str,
    )
    ap.add_argument(
        "-mu",
        "--bp_merge_units",
        help="Merge HOR units into HOR blocks within this number of base pairs.",
        type=int,
        default=DEF_BP_MERGE_UNITS,
    )
    ap.add_argument(
        "-mb",
        "--bp_merge_blks",
        help="Merge HOR blocks into HOR arrays within this number of bases pairs.",
        type=int,
        default=DEF_BP_MERGE_BLKS,
    )
    ap.add_argument(
        "-ub",
        "--min_blk_hor_units",
        help="HOR blocks must have at least n HOR units unbroken.",
        type=int,
        default=DEF_MIN_BLK_HOR_UNITS,
    )
    ap.add_argument(
        "-ua",
        "--min_arr_hor_units",
        help="Require that an HOR array have at least n HOR units.",
        type=int,
        default=DEF_MIN_ARR_HOR_UNITS,
    )
    ap.add_argument(
        "-fp",
        "--min_arr_prop",
        help="Require that an HOR array has at least this proportion of HOR by length.",
        type=float,
        default=DEF_MIN_ARR_PROP,
    )
    ap.add_argument(
        "-fl",
        "--min_arr_len",
        help="Require that an HOR array is this size in bp.",
        type=int,
        default=DEF_MIN_ARR_LEN,
    )
    return None


def calculate_hor_length(
    infile: TextIO,
    bp_merge_units: int,
    bp_merge_blks: int,
    min_blk_hor_units: int,
    min_arr_hor_units: int,
    min_arr_len: int,
    min_arr_prop: int,
    output: TextIO,
    output_strand: str | None = None,
    rmfile: TextIO | None = None,
) -> int:
    """
    Calculate HOR array length from HumAS-HMMER structural variation row output.

    ### Parameters
    `infile`
        Input bed file made from HumAS-HMMER output.
        Expects the following columns: `{chrom, chrom_st, chrom_end, hor, 0, strand, ...}`.
    `rmfile`
        Input RepeatMasker file.
        Used to prevent merging across other repeat types.
    `bp_merge_units`
        Merge HOR units into HOR blocks within this number of base pairs.
    `bp_merge_blks`
        Merge HOR blocks into HOR arrays within this number of bases pairs.
    `min_blk_hor_units`
        Grouped stv rows must have at least `n` HOR units unbroken.
    `min_arr_hor_units`
        Require that an HOR array have at least `n` HOR units.
    `min_arr_len`
        Require that an HOR array is this size in bp.
    `min_arr_prop`
        Require that an HOR array has at least this proportion of HORs by length.
    `output`
        Output bed file with HOR array lengths.
        Columns: `{chrom, chrom_st, chrom_end, length}`.
    `output_strand`
        Output bed file with HOR array lengths by strand.
        Columns: `{chrom, chrom_st, chrom_end, length, strand}`.

    ### Returns
    0 if successful.
    """
    df_stv = read_stv(infile)
    df_rm = read_rm(rmfile) if rmfile else None

    df_all_len, df_all_strand_len = hor_array_length(
        df_stv=df_stv,
        df_rm=df_rm,
        bp_merge_units=bp_merge_units,
        bp_merge_blks=bp_merge_blks,
        min_blk_hor_units=min_blk_hor_units,
        min_arr_hor_units=min_arr_hor_units,
        min_arr_len=min_arr_len,
        min_arr_prop=min_arr_prop,
        output_strand=isinstance(output_strand, str),
    )

    if output_strand:
        format_and_output_lengths(
            df_all_strand_len, output_strand, DEF_OUTPUT_BED_COLS_STRAND
        )

    format_and_output_lengths(df_all_len, output, DEF_OUTPUT_BED_COLS)
    return 0
