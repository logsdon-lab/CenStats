import re
import polars as pl

from .orientation import Orientation
from .constants import RGX_CHR


def get_contig_similarity_by_edit_dst(
    contigs: list[str],
    ref_contigs: list[str],
    edit_dst: list[int],
    orientation: list[Orientation],
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
    edit_distance_highest_dst_filter = pl.col("dst_perc") == pl.col("dst_perc").min()

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
                edit_distance_highest_dst_filter
            )

        if df_filter_ort_res_same_chr_grp.is_empty():
            df_filter_ort_res_same_chr_grp = df_edit_distance_res_same_chr_grp.filter(
                edit_distance_highest_dst_filter
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
