import re
import polars as pl
from typing import NamedTuple, Generator
from .constants import RGX_CHR


class RefCenContigs(NamedTuple):
    chr: str
    ref: str
    df: pl.DataFrame


def split_ref_rm_input_by_contig(
    df_ref: pl.DataFrame,
) -> Generator[tuple[str, RefCenContigs], None, None]:
    for ref, df_ref_grp in df_ref.group_by(["contig"]):
        ref = ref[0]
        mtch_ref_chr_name = re.search(RGX_CHR, ref)
        if not mtch_ref_chr_name:
            continue

        ref_chr_name = mtch_ref_chr_name.group()

        yield ref, RefCenContigs(ref_chr_name, ref, df_ref_grp)
