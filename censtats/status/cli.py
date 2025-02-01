import re
import sys
import argparse
import polars as pl
import editdistance
from loguru import logger
from typing import TextIO, TYPE_CHECKING, Any

from .repeat_jaccard_index import jaccard_index, get_contig_similarity_by_jaccard_index
from .repeat_edit_dst import get_contig_similarity_by_edit_dst
from .constants import (
    CHROMOSOMES_13_21,
    CHROMOSOMES_14_22,
    MAX_ALR_LEN_THR,
    RGX_CHR,
    EDGE_LEN,
    EDGE_PERC_ALR_THR,
    DST_PERC_THR,
)
from .orientation import Orientation
from .reference import split_ref_rm_input_by_contig
from .reader import read_repeatmasker_output
from .partial_cen import is_partial_centromere


if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def join_summarize_results(
    df_partial_contig_res: pl.DataFrame,
    df_jaccard_index_res: pl.DataFrame,
    df_edit_distance_res: pl.DataFrame,
    df_edit_distance_same_chr_res: pl.DataFrame,
    *,
    reference_prefix: str,
) -> pl.DataFrame:
    df_joined = (
        df_partial_contig_res.join(
            df_jaccard_index_res.join(df_edit_distance_res, on="contig")
            .group_by("contig")
            .first(),
            on="contig",
            how="left",
        )
        # Add default ort per contig.
        .join(df_edit_distance_same_chr_res, on="contig", how="left")
        .with_columns(ctg_chr=pl.col("contig").str.extract(RGX_CHR.pattern))
    )

    return (
        df_joined.select(
            contig=pl.col("contig"),
            # Extract chromosome name.
            # Both results must concur.
            final_contig=pl.when(pl.col("ref") == pl.col("ref_right"))
            .then(pl.col("ref").str.extract(RGX_CHR.pattern))
            .otherwise(pl.col("contig").str.extract(RGX_CHR.pattern)),
            # Only use orientation if both agree. Otherwise, replace with best same chr ort.
            reorient=pl.when(pl.col("ref") == pl.col("ref_right"))
            .then(pl.col("ort"))
            .otherwise(None)
            .fill_null(pl.col("ort_same_chr")),
            partial=pl.col("partial"),
        )
        # Replace chr name in original contig.
        .with_columns(
            final_contig=pl.col("contig").str.replace(
                RGX_CHR.pattern, pl.col("final_contig")
            ),
            # Never reorient if reference.
            reorient=pl.when(pl.col("contig").str.starts_with(reference_prefix))
            .then(
                pl.col("reorient").str.replace(
                    str(Orientation.Reverse), str(Orientation.Forward)
                )
            )
            .otherwise(pl.col("reorient")),
        )
        # Take only first row per contig.
        .group_by("contig", maintain_order=True)
        .first()
    )


def check_cens_status(
    input_rm: TextIO,
    output: TextIO,
    reference_rm: str,
    *,
    reference_prefix: str,
    dst_perc_thr: float = DST_PERC_THR,
    edge_perc_alr_thr: float = EDGE_PERC_ALR_THR,
    edge_len: int = EDGE_LEN,
    max_alr_len_thr: int = MAX_ALR_LEN_THR,
    restrict_13_21: bool = False,
    restrict_14_22: bool = False,
) -> int:
    df_ctg = read_repeatmasker_output(input_rm)
    df_ref = read_repeatmasker_output(reference_rm).filter(
        pl.col("contig").str.starts_with(reference_prefix)
    )

    contigs, refs, dsts, orts = [], [], [], []
    jcontigs, jrefs, jindex = [], [], []
    pcontigs, pstatus = [], []

    # Split ref dataframe by chromosome.
    df_ref_grps = dict(split_ref_rm_input_by_contig(df_ref))
    logger.info(f"Read {len(df_ref_grps)} reference dataframes.")

    for ctg, df_ctg_grp in df_ctg.group_by(["contig"]):
        ctg_name = ctg[0]
        logger.info(f"Evaluating {ctg_name} with {df_ctg_grp.shape[0]} repeats.")

        mtch_chr_name = re.search(RGX_CHR, ctg_name)
        if not mtch_chr_name:
            continue
        chr_name = mtch_chr_name.group()

        # Check if partial ctg.
        pcontigs.append(ctg_name)
        is_partial = is_partial_centromere(
            df_ctg_grp,
            edge_len=edge_len,
            edge_perc_alr_thr=edge_perc_alr_thr,
            max_alr_len_thr=max_alr_len_thr,
        )
        pstatus.append(is_partial)

        for ref_name, ref_ctg in df_ref_grps.items():
            df_ref_grp = ref_ctg.df

            # Special case for 13 and 21 and 14 and 22.
            if (
                (chr_name in CHROMOSOMES_13_21 and ref_ctg.chr not in CHROMOSOMES_13_21)
                and restrict_13_21
            ) or (
                (chr_name in CHROMOSOMES_14_22 and ref_ctg.chr not in CHROMOSOMES_14_22)
                and restrict_14_22
            ):
                continue

            dst_fwd = editdistance.eval(
                df_ref_grp["type"].to_list(),
                df_ctg_grp["type"].to_list(),
            )
            dst_rev = editdistance.eval(
                df_ref_grp["type"].to_list(),
                df_ctg_grp["type"].reverse().to_list(),
            )

            repeat_type_jindex = jaccard_index(
                set(df_ref_grp["type"]), set(df_ctg_grp["type"])
            )
            jcontigs.append(ctg_name)
            jrefs.append(ref_name)
            jindex.append(repeat_type_jindex)

            contigs.append(ctg_name)
            contigs.append(ctg_name)
            refs.append(ref_name)
            refs.append(ref_name)
            orts.append(str(Orientation.Forward))
            orts.append(str(Orientation.Reverse))
            dsts.append(dst_fwd)
            dsts.append(dst_rev)

    df_jaccard_index_res = get_contig_similarity_by_jaccard_index(
        jcontigs, jrefs, jindex
    )
    (
        df_filter_edit_distance_res,
        df_filter_ort_same_chr_res,
    ) = get_contig_similarity_by_edit_dst(
        contigs, refs, dsts, orts, dst_perc_thr=dst_perc_thr
    )
    df_partial_contig_res = pl.DataFrame({"contig": pcontigs, "partial": pstatus})

    res = join_summarize_results(
        df_partial_contig_res=df_partial_contig_res,
        df_jaccard_index_res=df_jaccard_index_res,
        df_edit_distance_res=df_filter_edit_distance_res,
        df_edit_distance_same_chr_res=df_filter_ort_same_chr_res,
        reference_prefix=reference_prefix,
    )

    res.write_csv(output, include_header=False, separator="\t")
    logger.info("Finished checking centromeres.")

    return 0


def add_status_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "status",
        description="Determines if centromeres are incorrectly oriented/mapped with respect to a reference.",
    )
    ap.add_argument(
        "-i",
        "--input",
        help="Input RepeatMasker output. Should contain contig reference. Expects no header.",
        type=argparse.FileType("rb"),
        required=True,
    )
    ap.add_argument(
        "-o",
        "--output",
        help="List of contigs with actions required to fix.",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )
    ap.add_argument(
        "-r",
        "--reference",
        required=True,
        type=str,
        help="Reference RM dataframe.",
    )
    ap.add_argument(
        "--dst_perc_thr",
        default=DST_PERC_THR,
        type=float,
        help="Edit distance percentile threshold. Lower is more stringent.",
    )
    ap.add_argument(
        "--edge_perc_alr_thr",
        default=EDGE_PERC_ALR_THR,
        type=float,
        help="Percent ALR on edges of contig to be considered a partial centromere.",
    )
    ap.add_argument(
        "--edge_len",
        default=EDGE_LEN,
        type=int,
        help="Edge len to calculate edge_perc_alr_thr.",
    )
    ap.add_argument(
        "--max_alr_len_thr",
        default=MAX_ALR_LEN_THR,
        type=int,
        help="Length of largest ALR needed in a contig to not be considered a partial centromere.",
    )
    ap.add_argument(
        "--reference_prefix", default="chm13", type=str, help="Reference prefix."
    )
    ap.add_argument(
        "--restrict_13_21",
        action="store_true",
        help="Restrict mapping to chromosomes 13 and 21 for chr13 and chr21 contigs.",
    )
    ap.add_argument(
        "--restrict_14_22",
        action="store_true",
        help="Restrict mapping to chromosomes 14 and 22 for chr14 and chr22 contigs.",
    )

    return None
