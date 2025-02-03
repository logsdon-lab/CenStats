from collections import Counter
import polars as pl
import intervaltree as it

from typing import Any, Callable, NamedTuple

from .constants import (
    DEF_BP_MERGE_BLKS,
    DEF_BP_MERGE_UNITS,
    DEF_MERGE_RBLACKLIST,
    DEF_MIN_ARR_HOR_UNITS,
    DEF_MIN_BLK_HOR_UNITS,
    DEF_MIN_ARR_PROP,
    DEF_MIN_ARR_LEN,
    DEF_OUTPUT_BED_COLS,
    DEF_OUTPUT_BED_COLS_STRAND,
)
from ..common import merge_itvs


class MergeHORData(NamedTuple):
    name: str
    hor_count: int
    hor_len: int


def bed_to_itree(df: pl.DataFrame, fn_mdata: Callable[[int, int], Any] | None = None):
    return it.IntervalTree(
        it.Interval(st, end, fn_mdata(st, end) if fn_mdata else None)
        for st, end in df.select("chrom_st", "chrom_end").iter_rows()
    )


def merge_hor_unit_itvs(i1: it.Interval, i2: it.Interval) -> it.Interval:
    i1_data: MergeHORData = i1.data
    i2_data: MergeHORData = i2.data
    # Keep count of number of intervals merged and merged length.
    new_data = MergeHORData(
        i1_data.name, i1_data.hor_count + 1, i1_data.hor_len + i2_data.hor_len
    )
    return it.Interval(i1.begin, i2.end, new_data)


def group_by_dst(df: pl.DataFrame, dst: int, group_name: str) -> pl.DataFrame:
    try:
        df = df.drop("index")
    except pl.exceptions.ColumnNotFoundError:
        pass
    return (
        df.with_columns(
            # c1  st1 (end1)
            # c1 (st2) end2
            dst_behind=(pl.col("chrom_st") - pl.col("chrom_end").shift(1)).fill_null(0),
            dst_ahead=(pl.col("chrom_st").shift(-1) - pl.col("chrom_end")).fill_null(0),
        )
        .with_row_index()
        .with_columns(
            **{
                # Group HOR units based on distance.
                group_name: pl.when(pl.col("dst_behind").le(dst))
                # We assign 0 if within merge dst.
                .then(pl.lit(0))
                # Otherwise, give unique index.
                .otherwise(pl.col("index") + 1)
                # Then create run-length ID to group on.
                # Contiguous rows within distance will be grouped together.
                .rle_id()
            },
        )
        .with_columns(
            # Adjust groups in scenarios where should join group ahead or behind but given unique group.
            # B:64617 A:52416 G:1
            # B:52416 A:1357  G:2 <- This should be group 3.
            # B:1357  A:1358  G:3
            pl.when(pl.col("dst_behind").le(dst) & pl.col("dst_ahead").le(dst))
            .then(pl.col(group_name))
            .when(pl.col("dst_behind").le(dst))
            .then(pl.col(group_name).shift(1))
            .when(pl.col("dst_ahead").le(dst))
            .then(pl.col(group_name).shift(-1))
            .otherwise(pl.col(group_name))
        )
    )


def hor_array_length(
    df_stv: pl.DataFrame,
    df_rm: pl.DataFrame | None = None,
    bp_merge_units: int = DEF_BP_MERGE_UNITS,
    bp_merge_blks: int = DEF_BP_MERGE_BLKS,
    min_blk_hor_units: int = DEF_MIN_BLK_HOR_UNITS,
    min_arr_hor_units: int = DEF_MIN_ARR_HOR_UNITS,
    min_arr_len: int = DEF_MIN_ARR_LEN,
    min_arr_prop: float = DEF_MIN_ARR_PROP,
    *,
    output_strand: bool = True,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    dfs: list[pl.DataFrame] = []
    dfs_strand: list[pl.DataFrame] = []
    for ctg_name, df_chr in df_stv.sort(by=["chrom", "chrom_st"]).group_by(
        ["chrom"], maintain_order=True
    ):
        ctg_name = ctg_name[0]
        df_live_hor = group_by_dst(
            df_chr.filter(pl.col("name").str.contains("L")),
            bp_merge_units,
            "live_group",
        ).filter(
            # Filter any live group with fewer than required number of HOR units.
            pl.col("live_group").count().over("live_group") >= min_blk_hor_units
        )
        itvs_live_hor = bed_to_itree(
            df_live_hor, lambda st, end: MergeHORData(ctg_name, 1, end - st)
        )

        # Filter to chrom
        if isinstance(df_rm, pl.DataFrame):
            itvs_chr_rm = it.IntervalTree(
                it.Interval(st, end, rtype)
                for st, end, rtype in df_rm.filter(pl.col("contig") == ctg_name)
                .select("start", "end", "type")
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
            return most_common_repeat not in DEF_MERGE_RBLACKLIST

        merged_itvs = merge_itvs(
            itvs_live_hor.iter(),
            dst=bp_merge_blks,
            fn_cmp=check_correct_merge,
            fn_merge_itv=merge_hor_unit_itvs,
        )
        if output_strand:
            # Group HOR units by strand into blocks
            # Requiring at least the min_blk_hor_units per strand block.
            df_live_hor = (
                group_by_dst(
                    df_live_hor.with_columns(
                        strand_group=pl.col("strand").rle_id()
                    ).filter(
                        pl.col("strand_group").count().over("strand_group")
                        >= min_blk_hor_units
                    ),
                    # Allow grouping by merge_blks as this is final group check before splitting by strand_group
                    bp_merge_blks,
                    "live_group",
                )
                # Take both strand and distance into consideration.
                .with_columns(
                    strand_group=pl.col("strand").rle_id() + pl.col("live_group")
                )
            )
            for _, df_strand_group in df_live_hor.group_by(["strand_group"]):
                strand = df_strand_group.get_column("strand")[0]
                itvs_strand = bed_to_itree(
                    df_strand_group, lambda st, end: MergeHORData(strand, 1, end - st)
                )
                merged_strand_itvs = merge_itvs(
                    itvs_strand.iter(),
                    dst=bp_merge_blks,
                    fn_cmp=check_correct_merge,
                    fn_merge_itv=merge_hor_unit_itvs,
                )
                df_strand = pl.DataFrame(
                    (
                        (
                            ctg_name,
                            itv.begin,
                            itv.end,
                            itv.length(),
                            itv.data[1],
                            itv.data[2] / itv.length(),
                            strand,
                        )
                        for itv in merged_strand_itvs
                    ),
                    orient="row",
                    schema=DEF_OUTPUT_BED_COLS_STRAND,
                ).filter(
                    (pl.col("score") >= min_arr_hor_units)
                    & (pl.col("prop") >= min_arr_prop)
                    & (pl.col("name") >= min_arr_len)
                )
                dfs_strand.append(df_strand)

        df = (
            pl.DataFrame(
                (
                    (
                        ctg_name,
                        itv.begin,
                        itv.end,
                        itv.length(),
                        itv.data[1],
                        itv.data[2] / itv.length(),
                    )
                    for itv in merged_itvs
                ),
                orient="row",
                schema=DEF_OUTPUT_BED_COLS,
            )
            # Require that array has at least n merged HOR units.
            .filter(
                (pl.col("score") >= min_arr_hor_units)
                & (pl.col("prop") >= min_arr_prop)
                & (pl.col("name") >= min_arr_len)
            )
        )

        dfs.append(df)

    df_all = pl.concat(dfs).sort(by=["chrom", "chrom_st"])
    df_all_strand = (
        pl.concat(dfs_strand)
        if output_strand
        else pl.DataFrame(schema=DEF_OUTPUT_BED_COLS_STRAND)
    ).sort(by=["chrom", "chrom_st"])
    return df_all, df_all_strand
