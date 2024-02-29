import polars as pl

from .orientation import Orientation


def determine_acro_p_arm_ort(df: pl.DataFrame) -> tuple[Orientation, pl.DataFrame]:
    """
    Map the p-arm of an acrocentric chromosome centromeric contig to an orientation.

    ### Args
    `df`
        RepeatMasker annotation dataframe of a single centromeric contig.

    ### Returns
    `Orientation` and a `pl.DataFrame` of repeats from the largest ALR.
    """
    start_bp_pos = df.row(0, named=True)["end"]
    end_bp_pos = df.row(-1, named=True)["end"]
    largest_alr_repeat = df.filter(
        (pl.col("dst") == pl.col("dst").max().over(pl.col("type")))
        & (pl.col("type") == "ALR/Alpha")
    )
    largest_alr_repeat_mdpt_pos = (
        largest_alr_repeat["end"][0] - largest_alr_repeat["start"][0]
    )
    abs_dst_to_start = abs(start_bp_pos - largest_alr_repeat_mdpt_pos)
    abs_dst_to_end = abs(end_bp_pos - largest_alr_repeat_mdpt_pos)

    if abs_dst_to_start < abs_dst_to_end:
        return (Orientation.Forward, largest_alr_repeat)
    else:
        return (Orientation.Reverse, largest_alr_repeat)


def get_q_arm_acro_chr(df_ctg_grp: pl.DataFrame) -> pl.DataFrame:
    """
    Get the q-arm of an acrocentric chromosome's centromere.

    ### Args
    `df`
        RepeatMasker annotation dataframe of a single centromeric contig.

    ### Returns
    `pl.DataFrame` of repeats of the q-arm of acrocentric chromosome.
    """
    p_arm_ort, alr_repeat = determine_acro_p_arm_ort(df_ctg_grp)
    # | p | alr | q |
    if p_arm_ort == Orientation.Forward:
        return df_ctg_grp.filter(pl.col("start") > alr_repeat["end"][0])
    # | q | alr | p |
    else:
        return df_ctg_grp.filter(pl.col("end") > alr_repeat["start"][0])


def flatten_repeats(df: pl.DataFrame, *, window_size=5) -> pl.DataFrame:
    """
    Flattens/denoises sequences of repeats by:
    1. Creating overlapping repeat windows of `window_size`.
    2. Finding the largest repeat in that window and setting the current row type to that type.
    3. Grouping sequences of repeats of the same type.
    4. Merging the group and recalculating the start, end, and distance of the merged repeat.

    This helps find the position of the HOR array.

    ### Args
    `df`
        RepeatMasker annotation dataframe of a single centromeric contig.

    ### Returns
    Flattened `pl.DataFrame` of repeats with the columns:
    1. `start`
    2. `end`
    3. `type`
    4. `dst`
    """
    return (
        df.with_columns(
            type=pl.Series(
                (
                    df.slice(i, window_size)
                    .select(pl.col("type").filter(pl.col("dst") == pl.col("dst").max()))
                    .get_column("type")[0]
                )
                for i in range(len(df["type"]))
            )
        )
        .with_columns(id=pl.col("type").rle_id())
        .group_by("id")
        .agg(
            [
                pl.col("start").min().alias("start"),
                pl.col("end").max().alias("end"),
                pl.col("type").first(),
            ]
        )
        .drop("id")
        .with_columns(dst=pl.col("end") - pl.col("start"))
        .sort(by="start")
    )
