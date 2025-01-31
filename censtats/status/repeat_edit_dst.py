import re
import polars as pl

from typing import Iterator, Optional

from .constants import RGX_CHR, REPEAT_SPLIT_LEN


def split_repeats(x: int, div: int) -> Iterator[int]:
    """
    Explodes/expands repeats per div.
    * ex. div = 1000
            * 2001 bp ALR = [1000 bp ALR, 1000 bp ALR, 1 bp ALR]
    """
    d, m = divmod(x, div)
    for d in range(d):
        yield div
    yield m


def expand_repeat_dst(
    df_ctg_grp: pl.DataFrame,
    *,
    repeat_filter: Optional[pl.Expr] = None,
    bp_repeat_split: int = REPEAT_SPLIT_LEN,
) -> pl.DataFrame:
    """
    Expand the repeat distance.

    ### Args
    `df_ctg_grp`
        RepeatMasker annotation dataframe of a single centromeric contig.
    `repeat_filter`
        Expression to filter and expand a subset of repeats.
    `bp_repeat_split`
        Number of base pairs to split the repeat by.

    ### Returns
    `pl.DataFrame` of expanded repeats.
    """
    return df_ctg_grp.with_columns(
        pl.when(repeat_filter if repeat_filter is not None else False)
        .then(
            pl.col("dst").map_elements(
                lambda x: list(split_repeats(x, bp_repeat_split)),
                return_dtype=pl.List(pl.Int64),
            )
        )
        .otherwise(pl.col("dst").cast(pl.List(pl.Int64)))
    ).explode("dst")


def get_contig_similarity_by_edit_dst(
    contigs: list[str],
    ref_contigs: list[str],
    edit_dst: list[int],
    orientation: list[str],
    *,
    dst_perc_thr: float,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    df_edit_distance_res = pl.DataFrame(
        {"contig": contigs, "ref": ref_contigs, "dst": edit_dst, "ort": orientation}
    ).with_columns(
        dst_perc=(pl.col("dst").rank() / pl.col("dst").count()).over("contig"),
    )

    # Filter results so only:
    # * Matches gt x percentile.
    # * Distances lt y percentile.
    dfs_filtered_edit_distance_res = []
    dfs_filtered_ort_same_chr_res = []

    edit_distance_thr_filter = pl.col("dst_perc") < dst_perc_thr
    edit_distance_lowest_dst_filter = pl.col("dst_perc") == pl.col("dst_perc").min()

    for contig, df_edit_distance_res_grp in df_edit_distance_res.group_by(["contig"]):
        mtch_chr_name = re.search(RGX_CHR, contig[0])
        if not mtch_chr_name:
            continue
        chr_name = mtch_chr_name.group()

        # Only look at same chr to determine default ort.
        df_edit_distance_res_same_chr_grp = df_edit_distance_res_grp.filter(
            pl.col("ref").str.contains(f"{chr_name}:")
        )

        df_filter_edit_distance_res_grp = df_edit_distance_res_grp.filter(
            edit_distance_thr_filter
        )
        df_filter_ort_res_same_chr_grp = df_edit_distance_res_same_chr_grp.filter(
            edit_distance_thr_filter
        )

        # If none found, default to highest number of matches.
        if df_filter_edit_distance_res_grp.is_empty():
            df_filter_edit_distance_res_grp = df_edit_distance_res_grp.filter(
                edit_distance_lowest_dst_filter
            )

        if df_filter_ort_res_same_chr_grp.is_empty():
            df_filter_ort_res_same_chr_grp = df_edit_distance_res_same_chr_grp.filter(
                edit_distance_lowest_dst_filter
            )

        dfs_filtered_edit_distance_res.append(df_filter_edit_distance_res_grp)
        dfs_filtered_ort_same_chr_res.append(df_filter_ort_res_same_chr_grp)

    df_filter_edit_distance_res: pl.DataFrame = pl.concat(
        dfs_filtered_edit_distance_res
    )
    df_filter_ort_same_chr_res: pl.DataFrame = pl.concat(dfs_filtered_ort_same_chr_res)

    df_filter_edit_distance_res = (
        df_filter_edit_distance_res
        # https://stackoverflow.com/a/74336952
        .with_columns(pl.col("dst").min().over("contig").alias("lowest_dst"))
        .filter(pl.col("dst") == pl.col("lowest_dst"))
        .select(["contig", "ref", "dst", "ort"])
    )
    # Get pair with lowest dst to get default ort.
    df_filter_ort_same_chr_res = (
        df_filter_ort_same_chr_res.with_columns(
            pl.col("dst").min().over("contig").alias("lowest_dst")
        )
        .filter(pl.col("dst") == pl.col("lowest_dst"))
        .select(["contig", "ort"])
        .rename({"ort": "ort_same_chr"})
    )

    return df_filter_edit_distance_res, df_filter_ort_same_chr_res
