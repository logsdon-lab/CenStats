import polars as pl

from .orientation import Orientation


def determine_acro_arm_ort(df: pl.DataFrame) -> tuple[Orientation, pl.DataFrame]:
    """
    Map the arms of an acrocentric chromosome centromeric contig to an orientation.

    ### Args
    `df`
        RepeatMasker annotation dataframe of a single centromeric contig.

    ### Returns
    `Orientation` and the row of the ALR.
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


def get_p_arm_acro_chr(df_ctg_grp: pl.DataFrame, *, additional_bp=500_000):
    arm_ort, alr_repeat = determine_acro_arm_ort(df_ctg_grp)
    if arm_ort == Orientation.Forward:
        return df_ctg_grp.filter(pl.col("end") < alr_repeat["end"][0] + additional_bp)
    else:
        return df_ctg_grp.filter(
            pl.col("start") > alr_repeat["start"][0] - additional_bp
        )


def flatten_repeats(df: pl.DataFrame) -> pl.DataFrame:
    """
    Flattens/denoises sequences of repeats by:
    1. Removing single repeats sandwiched in-between two repeats of the same type.
    2. Grouping sequences of repeats of the same type.
    3. Merging the group and recalculating the start, end, and distance of the merged repeat.

    This helps find the position of the HOR array.

    ### Args
    `df`
        RepeatMasker annotation dataframe of a single centromeric contig.

    ### Returns
    Flattened dataframe of repeats with the columns:
    1. `start`
    2. `end`
    3. `type`
    4. `dst`
    """
    return (
        df.with_columns(
            type=pl.when(pl.col("type").shift(n=-1) == pl.col("type").shift(n=1))
            .then(pl.col("type").shift(n=1))
            .otherwise(pl.col("type")),
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
