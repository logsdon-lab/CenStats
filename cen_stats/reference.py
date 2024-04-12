import re
import polars as pl
from typing import NamedTuple, Generator
from .acrocentrics import flatten_repeats, get_q_arm_acro_chr
from .constants import ACROCENTRIC_CHROMOSOMES, HOR_LEN_THR, RGX_CHR


class RefCenContigs(NamedTuple):
    ref: str
    df: pl.DataFrame
    flat_df: pl.DataFrame
    num_hor_arrays: int


def split_ref_rm_input_by_contig(
    df_ref: pl.DataFrame,
) -> Generator[tuple[str, RefCenContigs], None, None]:
    for ref, df_ref_grp in df_ref.group_by(["contig"]):
        ref = ref[0]
        mtch_ref_chr_name = re.search(RGX_CHR, ref)
        if not mtch_ref_chr_name:
            continue

        ref_chr_name = mtch_ref_chr_name.group()

        df_ref_flatten_grp = flatten_repeats(df_ref_grp)
        num_hor_arrays = len(
            df_ref_flatten_grp.filter(
                (pl.col("type") == "ALR/Alpha") & (pl.col("dst") > HOR_LEN_THR)
            )
        )
        # Also adjust for reference acrocentrics.
        if ref_chr_name in ACROCENTRIC_CHROMOSOMES:
            df_ref_flatten_grp = get_q_arm_acro_chr(df_ref_flatten_grp)

        yield ref, RefCenContigs(ref, df_ref_grp, df_ref_flatten_grp, num_hor_arrays)