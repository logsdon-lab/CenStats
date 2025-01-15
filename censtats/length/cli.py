import re
import sys
import argparse
import polars as pl

from typing import TYPE_CHECKING, Any, TextIO

from .constants import (
    RGX_CHR,
    DEF_BP_JUMP_LEN_THR,
    DEF_ARR_LEN_THR,
    HOR_BP_LEN,
    DEF_EXP_STV_ROW_BED_COLS,
    DEF_OUTPUT_BED_COLS,
)


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
        "--input",
        help=f"Input stv row bed file produced by HumAS-HMMER and stv. Expects columns: {DEF_EXP_STV_ROW_BED_COLS}",
        type=argparse.FileType("rb"),
    )
    ap.add_argument(
        "-o",
        "--output",
        help=f"Output bed file with columns: {DEF_OUTPUT_BED_COLS}",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )
    ap.add_argument(
        "--bp_jump_thr",
        help="Base pair jump threshold to group by",
        type=int,
        default=DEF_BP_JUMP_LEN_THR,
    )
    ap.add_argument(
        "--arr_len_thr",
        help="Length threshold to filter out.",
        type=int,
        default=DEF_ARR_LEN_THR,
    )
    return None


def calculate_hor_length(
    infile: TextIO, bp_jump_thr: int, arr_len_thr: int, output: str
) -> int:
    """
    Calculate HOR array length from HumAS-HMMER structural variation row output.

    ### Parameters
    `infile`
        Input bed file made from HumAS-HMMER output.
        Expects the following columns: `{chr, start, stop, hor, 0, strand, ...}`.
    `bp_jump_thr`
        Base pair jump threshold to group by.
    `arr_len_thr`
        Length threshold of HOR array to filter out.
    `output`
        Output bed file with HOR array lengths.
        Columns: `{chr_name, start_pos, stop_pos, len}`.

    ### Returns
    0 if successful.
    """
    df = pl.read_csv(
        infile,
        separator="\t",
        columns=[0, 1, 2, 3, 4, 5],
        new_columns=DEF_EXP_STV_ROW_BED_COLS,
        has_header=False,
    )

    dfs = []
    for ctg_name, df_chr in df.group_by(["chr"], maintain_order=True):
        df_chr = df_chr.with_columns(len=pl.col("stop") - pl.col("start")).with_columns(
            mer=(pl.col("len") / HOR_BP_LEN).round()
        )
        ctg_name = ctg_name[0]
        mtch_chr_name = re.search(RGX_CHR, ctg_name)
        if mtch_chr_name is None:
            continue

        chr_name = mtch_chr_name.group(0)
        df_live_hor = df_chr.filter(pl.col("hor").str.contains("L"))

        # Specific edge case for chr8.
        if chr_name == "chr8" or chr_name == "chr10" or chr_name == "chr16":
            bp_jump_thr = 10_000
        elif chr_name == "chrY":
            bp_jump_thr = 2_000
        else:
            bp_jump_thr = bp_jump_thr

        df_bp_jumps = df_live_hor.with_columns(
            diff=pl.col("start") - pl.col("stop").shift(1)
        ).filter(pl.col("diff") > bp_jump_thr)

        if df_bp_jumps.is_empty():
            adj_start = df_live_hor.get_column("start").min()
            adj_stop = df_live_hor.get_column("stop").max()
            adj_len = adj_stop - adj_start

            if adj_len < arr_len_thr:
                continue

            dfs.append(
                pl.DataFrame(
                    {
                        "chr_name": ctg_name,
                        "start_pos": adj_start,
                        "stop_pos": adj_stop,
                        "len": adj_len,
                    }
                )
            )
            continue

        starts, stops = [], []
        for i, row in enumerate(df_bp_jumps.iter_rows()):
            prev_row = pl.DataFrame() if i == 0 else df_bp_jumps.slice(i - 1)
            next_row = df_bp_jumps.slice(i + 1)

            if prev_row.is_empty():
                starts.append(df_chr.get_column("start").min())
                stops.append(
                    df_chr.filter(pl.col("start") < row[1]).row(-1, named=True)["stop"]
                )

            if next_row.is_empty():
                starts.append(row[1])
                stops.append(df_chr.get_column("stop").max())
            else:
                starts.append(row[1])
                stops.append(
                    df_chr.filter(
                        pl.col("start") < next_row.get_column("start")[0]
                    ).row(-1, named=True)["stop"]
                )

        lens = []
        chr_mer_filter = None
        if chr_name == "chr10" or chr_name == "chr20":
            chr_mer_filter = pl.col("mer") >= 5
        elif chr_name == "chrY":
            chr_mer_filter = pl.col("mer") >= 30
        elif chr_name == "chr17":
            chr_mer_filter = pl.col("mer") >= 4

        for start, stop in zip(starts, stops):
            df_slice = (
                df_chr.filter(pl.col("start") >= start, pl.col("stop") <= stop)
                .with_columns(bp_jump=pl.col("start") - pl.col("stop").shift(1))
                .fill_null(0)
            )
            # Filter out mers based on chr.
            if chr_mer_filter is not None:
                df_slice = df_slice.filter(chr_mer_filter)

            if df_slice.is_empty():
                lens.append(0)
                continue
            df_slice_dst = (
                # df_slice.with_columns(len=pl.col("stop") - pl.col("start")).get_column("len").sum()
                df_slice.get_column("stop").max() - df_slice.get_column("start").min()
            )
            lens.append(df_slice_dst)

        lf = pl.LazyFrame(
            {
                "chr_name": ctg_name,
                "start_pos": starts,
                "stop_pos": stops,
                "len": lens,
            }
        )
        if (
            chr_name == "chr8"
            or chr_name == "chr10"
            or chr_name == "chr17"
            or chr_name == "chrY"
        ):
            arr_len_thr = 100_000
        else:
            arr_len_thr = arr_len_thr

        dfs.append(lf.filter(pl.col("len") > arr_len_thr).collect())

    df_all_dsts: pl.DataFrame = pl.concat(dfs)
    (
        df_all_dsts.with_columns(
            sort_idx=pl.col("chr_name")
            .str.extract("chr([0-9XY]+)")
            .replace({"X": "23", "Y": "24"})
            .cast(pl.Int32)
        )
        .sort(by="sort_idx")
        .select(DEF_OUTPUT_BED_COLS)
        .write_csv(output, include_header=False, separator="\t")
    )
    return 0
