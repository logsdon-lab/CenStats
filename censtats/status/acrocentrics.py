import polars as pl

from typing import Iterator, NamedTuple

from .constants import REPEAT_SPLIT_LEN


class AcroArms(NamedTuple):
    p_arm: pl.DataFrame
    q_arm: pl.DataFrame


def determine_acro_arms(df: pl.DataFrame, *, required_num_rtypes: int = 5) -> AcroArms:
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
        largest_alr_repeat["dst"][0] / 2
    ) + largest_alr_repeat["start"][0]
    abs_dst_to_start = abs(start_bp_pos - largest_alr_repeat_mdpt_pos)
    abs_dst_to_end = abs(end_bp_pos - largest_alr_repeat_mdpt_pos)
    df_leftarm = df.filter(pl.col("end") < largest_alr_repeat["start"][0])
    df_rightarm = df.filter(pl.col("start") > largest_alr_repeat["end"][0])
    l_rtypes = df_leftarm["type"].unique()
    r_rtypes = df_rightarm["type"].unique()

    # Assumption: If less than required_num_rtypes different repeat types, then it is a partial centromere with a break at the checked arm.
    if len(l_rtypes) < required_num_rtypes:
        return AcroArms(p_arm=df_leftarm, q_arm=df_rightarm)
    elif len(r_rtypes) < required_num_rtypes:
        return AcroArms(p_arm=df_rightarm, q_arm=df_leftarm)

    # Check for the orientation of the p-arm.
    # Assumption: p-arm is the larger arm of the two.
    # Required to have a minimum of required_num_rtypes different repeat types to be considered correct. Could be partial.
    if abs_dst_to_start > abs_dst_to_end:
        # | p | alr | q |
        df_parm = df_leftarm
        df_qarm = df_rightarm
    else:
        # | q | alr | p |
        df_parm = df_rightarm
        df_qarm = df_leftarm

    return AcroArms(df_parm, df_qarm)


def get_q_arm_acro_chr(
    df_ctg_grp: pl.DataFrame, *, bp_repeat_split: int = REPEAT_SPLIT_LEN
) -> pl.DataFrame:
    """
    Get the q-arm of an acrocentric chromosome's centromere.

    ### Args
    `df_ctg_grp`
        RepeatMasker annotation dataframe of a single centromeric contig.
    `bp_repeat_split`
        Number of base pairs to split the repeat by.

    ### Returns
    `pl.DataFrame` of repeats of the q-arm of acrocentric chromosome.
    """
    acro_arms = determine_acro_arms(df_ctg_grp)

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

    return acro_arms.q_arm.with_columns(
        pl.col("dst").map_elements(
            lambda x: list(split_repeats(x, bp_repeat_split)),
            return_dtype=pl.List(pl.Int64),
        )
    ).explode("dst")
