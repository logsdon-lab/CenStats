import sys
import argparse
import polars as pl
import intervaltree as it

from typing import TYPE_CHECKING, Any, TextIO
from collections import Counter

from .constants import (
    DEF_MIN_HOR_MONS,
    DEF_MIN_ARR_HOR_UNITS,
    DEF_MIN_GRP_HOR_UNITS,
    DEF_BP_MERGE,
    DEF_EXP_STV_ROW_BED_COLS,
    DEF_OUTPUT_BED_COLS,
    DEF_OUTPUT_BED_COLS_STRAND,
    DEF_RM_COLS,
    RM_COL_IDX,
    DEF_MERGE_RCLASSES,
    MON_LEN,
)
from ..common import merge_itvs


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
        help=f"Input stv row bed file produced by HumAS-HMMER and stv. Expects columns: {DEF_EXP_STV_ROW_BED_COLS}",
        type=argparse.FileType("rb"),
    )
    ap.add_argument(
        "-r",
        "--input_rm",
        help=f"Input RepeatMasker file with no header. Prevents joining across Expects columns: {DEF_RM_COLS}",
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
        "-m",
        "--bp_merge",
        help="Base pair merge threshold.",
        type=int,
        default=DEF_BP_MERGE,
    )
    ap.add_argument(
        "-g",
        "--min_grp_hor_units",
        help="Grouped stv rows must have at least n HOR units unbroken.",
        type=int,
        default=DEF_MIN_GRP_HOR_UNITS,
    )
    ap.add_argument(
        "-l",
        "--min_hor_mons",
        help="Require that each HOR have at least n monomers.",
        type=int,
        default=DEF_MIN_HOR_MONS,
    )
    ap.add_argument(
        "-a",
        "--min_arr_hor_units",
        help="Require that an HOR array have at least n HOR units.",
        type=int,
        default=DEF_MIN_ARR_HOR_UNITS,
    )
    return None


def bed_to_itree(df: pl.DataFrame, mdata: Any | None = None):
    return it.IntervalTree(
        it.Interval(st, end, mdata)
        for st, end in df.select("chrom_st", "chrom_end").iter_rows()
    )


def format_and_output_lengths(
    df: pl.DataFrame,
    output: TextIO | str,
    output_cols: list[str],
):
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


def calculate_hor_length(
    infile: TextIO,
    bp_merge: int,
    min_grp_hor_units: int,
    min_hor_mons: int,
    min_arr_hor_units: int,
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
    `bp_merge`
        Base pair merge threshold.
    `min_grp_hor_units`
        Grouped stv rows must have at least `n` HOR units unbroken.
    `min_hor_mons`
        Require that each HOR have at least `n` monomers.
    `min_arr_hor_units`
        Require that an HOR array have at least `n` HOR units.
    `output`
        Output bed file with HOR array lengths.
        Columns: `{chrom, chrom_st, chrom_end, length}`.
    `output_strand`
        Output bed file with HOR array lengths by strand.
        Columns: `{chrom, chrom_st, chrom_end, length, strand}`.

    ### Returns
    0 if successful.
    """
    df_stv = (
        pl.read_csv(
            infile,
            separator="\t",
            columns=[0, 1, 2, 3, 4, 5],
            new_columns=DEF_EXP_STV_ROW_BED_COLS,
            has_header=False,
        )
        .with_columns(
            mer=((pl.col("chrom_end") - pl.col("chrom_st")) / MON_LEN).round()
        )
        # Filter HOR units with fewer than required monomers.
        .filter(pl.col("mer") >= min_hor_mons)
    )
    if rmfile:
        df_rm = (
            pl.read_csv(
                rmfile,
                separator="\t",
                has_header=False,
                columns=RM_COL_IDX,
                new_columns=DEF_RM_COLS,
                truncate_ragged_lines=True,
            )
            .with_columns(
                ctg_name=pl.col("contig").str.extract(r"^(.*?):|^(.*?)$"),
                ctg_st=pl.col("contig")
                .str.extract(r":(\d+)-")
                .cast(pl.Int64)
                .fill_null(0),
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
    else:
        df_rm = None

    dfs = []
    dfs_strand = []
    for ctg_name, df_chr in df_stv.sort(by=["chrom", "chrom_st"]).group_by(
        ["chrom"], maintain_order=True
    ):
        ctg_name = ctg_name[0]
        df_live_hor = (
            df_chr
            # Only live regions.
            .filter(pl.col("name").str.contains("L"))
            .with_columns(
                len=pl.col("chrom_end") - pl.col("chrom_st"),
                # c1  st1 (end1)
                # c1 (st2) end2
                dst_behind=(
                    pl.col("chrom_st") - pl.col("chrom_end").shift(1)
                ).fill_null(0),
                dst_ahead=(
                    pl.col("chrom_st").shift(-1) - pl.col("chrom_end")
                ).fill_null(0),
            )
            .with_row_index()
            .with_columns(
                # Group HOR units based on distance.
                live_group=pl.when(pl.col("dst_behind").le(bp_merge))
                # We assign 0 if within merge dst.
                .then(pl.lit(0))
                # Otherwise, give unique index.
                .otherwise(pl.col("index") + 1)
                # Then create run-length ID to group on.
                # Contiguous rows within distance will be grouped together.
                .rle_id(),
            )
            .with_columns(
                # Adjust groups in scenarios where should join group ahead or behind but given unique group.
                # B:64617 A:52416 G:1
                # B:52416 A:1357  G:2 <- This should be group 3.
                # B:1357  A:1358  G:3
                pl.when(
                    pl.col("dst_behind").le(bp_merge) & pl.col("dst_ahead").le(bp_merge)
                )
                .then(pl.col("live_group"))
                .when(pl.col("dst_behind").le(bp_merge))
                .then(pl.col("live_group").shift(1))
                .when(pl.col("dst_ahead").le(bp_merge))
                .then(pl.col("live_group").shift(-1))
                .otherwise(pl.col("live_group"))
            )
            .filter(
                # Filter any live group with fewer than required number of HOR units.
                pl.col("live_group").count().over("live_group") >= min_grp_hor_units
            )
        )
        itvs_live_hor = bed_to_itree(df_live_hor, (ctg_name, 1))

        # Filter to chrom
        if isinstance(df_rm, pl.DataFrame):
            itvs_chr_rm = it.IntervalTree(
                it.Interval(st, end, rtype)
                for st, end, rtype in df_rm.filter(pl.col("contig") == ctg_name)
                .select("start", "end", "rClass")
                .iter_rows()
            )
        else:
            itvs_chr_rm = None

        # Then check what's between our intervals before merging.
        def check_correct_merge(itv_1: it.Interval, itv_2: it.Interval) -> bool:
            # If no repeatmasker tracks, no second check.
            if not itvs_chr_rm:
                return True
            itv_between = it.Interval(itv_1.end, itv_2.begin)
            # Prevent costly itree lookup if only small interval.
            if itv_between.length() <= 1:
                return True

            rm_ovl = itvs_chr_rm.overlap(itv_between)

            if not rm_ovl:
                return True
            repeat_count: Counter[str] = Counter()
            for ovl in rm_ovl:
                repeat_count[ovl.data] += ovl.overlap_size(itv_between)
            # Should never KeyError as we return early if rm_ovl is empty.
            most_common_repeat, _ = repeat_count.most_common(1)[0]

            # Only merge if most common repeat is allowed.
            return most_common_repeat in DEF_MERGE_RCLASSES

        merged_itvs = merge_itvs(
            itvs_live_hor.iter(),
            dst=bp_merge,
            fn_cmp=check_correct_merge,
            # Keep count of number of intervals merged.
            fn_merge_itv=lambda i1, i2: it.Interval(
                i1.begin, i2.end, (i1.data[0], i1.data[1] + 1)
            ),
        )
        if output_strand:
            df_live_hor = df_live_hor.with_columns(
                strand_group=pl.col("strand").rle_id()
            )
            for _, df_strand_group in df_live_hor.group_by(["strand_group"]):
                strand = df_strand_group.get_column("strand")[0]
                itvs_strand = bed_to_itree(df_strand_group, (strand, 1))
                merged_strand_itvs = merge_itvs(
                    itvs_strand.iter(),
                    dst=bp_merge,
                    fn_cmp=check_correct_merge,
                    fn_merge_itv=lambda i1, i2: it.Interval(
                        i1.begin, i2.end, (i1.data[0], i1.data[1] + 1)
                    ),
                )
                df_strand = pl.DataFrame(
                    (
                        (
                            ctg_name,
                            itv.begin,
                            itv.end,
                            itv.length(),
                            itv.data[1],
                            strand,
                        )
                        for itv in merged_strand_itvs
                    ),
                    orient="row",
                    schema=[
                        "chrom",
                        "chrom_st",
                        "chrom_end",
                        "name",
                        "score",
                        "strand",
                    ],
                ).filter(pl.col("score") >= min_arr_hor_units)
                dfs_strand.append(df_strand)

        df = (
            pl.DataFrame(
                ((ctg_name, itv.begin, itv.end, itv.data[1]) for itv in merged_itvs),
                orient="row",
                schema=["chrom", "chrom_st", "chrom_end", "score"],
            )
            .with_columns(name=pl.col("chrom_end") - pl.col("chrom_st"))
            # Require that array has at least n merged HOR units.
            .filter(pl.col("score") >= min_arr_hor_units)
        )

        dfs.append(df)

    if output_strand:
        df_all_strand_dsts = pl.DataFrame = pl.concat(dfs_strand)
        format_and_output_lengths(
            df_all_strand_dsts, output_strand, DEF_OUTPUT_BED_COLS_STRAND
        )

    df_all_dsts: pl.DataFrame = pl.concat(dfs)
    format_and_output_lengths(df_all_dsts, output, DEF_OUTPUT_BED_COLS)
    return 0
