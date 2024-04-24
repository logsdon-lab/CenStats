import re
import sys
import argparse
import polars as pl
import editdistance
from loguru import logger
from typing import TextIO

from .orientation import Orientation
from .repeat_jaccard_index import jaccard_index, get_contig_similarity_by_jaccard_index
from .repeat_edit_dst import get_contig_similarity_by_edit_dst
from .acrocentrics import get_q_arm_acro_chr, flatten_repeats
from .constants import (
    ACROCENTRIC_CHROMOSOMES,
    RGX_CHR,
    EDGE_LEN,
    EDGE_PERC_ALR_THR,
    DST_PERC_THR,
    HOR_LEN_THR,
)
from .reference import split_ref_rm_input_by_contig
from .reader import read_repeatmasker_output
from .partial_cen import is_partial_centromere


def join_summarize_results(
    df_partial_contig_res: pl.DataFrame,
    df_jaccard_index_res: pl.DataFrame,
    df_edit_distance_res: pl.DataFrame,
    df_edit_distance_same_chr_res: pl.DataFrame,
    *,
    reference_prefix: str,
) -> pl.DataFrame:
    return (
        # Join result dfs.
        # use partial contig res so always get all contigs.
        df_partial_contig_res.join(
            df_jaccard_index_res.join(df_edit_distance_res, on="contig")
            .group_by("contig")
            .first(),
            on="contig",
            how="left",
        )
        # Add default ort per contig.
        .join(df_edit_distance_same_chr_res, on="contig", how="left")
        .select(
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
                pl.col("reorient").str.replace(Orientation.Reverse, Orientation.Forward)
            )
            .otherwise(pl.col("reorient")),
        )
        # Take only first row per contig.
        .group_by("contig", maintain_order=True)
        .first()
    )


def check_cens_status(
    input_rm: str,
    output: TextIO,
    reference_rm: str,
    *,
    reference_prefix: str,
    dst_perc_thr: float = DST_PERC_THR,
    edge_perc_alr_thr: float = EDGE_PERC_ALR_THR,
    edge_len: int = EDGE_LEN,
    max_alr_len_thr: int = HOR_LEN_THR,
) -> int:
    df_ctg = read_repeatmasker_output(input_rm).collect()
    df_ref = (
        read_repeatmasker_output(reference_rm)
        .filter(pl.col("contig").str.starts_with(reference_prefix))
        .collect()
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
        pstatus.append(
            is_partial_centromere(
                df_ctg_grp,
                edge_len=edge_len,
                edge_perc_alr_thr=edge_perc_alr_thr,
                max_alr_len_thr=max_alr_len_thr,
            )
        )
        df_flatten_ctg_grp = flatten_repeats(df_ctg_grp)
        ctg_num_hor_arrays = len(
            df_flatten_ctg_grp.filter(
                (pl.col("type") == "ALR/Alpha") & (pl.col("dst") > HOR_LEN_THR)
            )
        )

        # For acros (13, 14, 15, 21, 21)
        # Adjust metrics to only use q-arm of chr.
        if chr_name in ACROCENTRIC_CHROMOSOMES:
            df_q_arm_ctg_grp = get_q_arm_acro_chr(df_flatten_ctg_grp)

        for ref_name, ref_ctg in df_ref_grps.items():
            # Check difference in number of HOR arrays between two contigs to determine if really acrocentric chr.
            # If diff in num of HOR arrays less than 3, assume same chr and align to q-arm.
            if (
                chr_name in ACROCENTRIC_CHROMOSOMES
                and ref_name in ACROCENTRIC_CHROMOSOMES
                and abs(ctg_num_hor_arrays - ref_ctg.num_hor_arrays) < 3
            ):
                # The flat_df is also just the q-arm
                df_ref_grp = ref_ctg.flat_df
                df_ctg_grp = df_q_arm_ctg_grp
            else:
                df_ref_grp = ref_ctg.df
                df_ctg_grp = df_ctg_grp

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
            orts.append(Orientation.Forward)
            orts.append(Orientation.Reverse)
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


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Determines if centromeres are incorrectly oriented/mapped with respect to a reference."
    )
    ap.add_argument(
        "-i",
        "--input",
        help="Input RepeatMasker output. Should contain contig reference. Expects no header.",
        type=str,
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
        default=HOR_LEN_THR,
        type=int,
        help="Length of largest ALR needed in a contig to not be considered a partial centromere.",
    )
    ap.add_argument(
        "--reference_prefix", default="chm13", type=str, help="Reference prefix."
    )
    args = ap.parse_args()

    return check_cens_status(
        args.input,
        args.output,
        args.reference,
        reference_prefix=args.reference_prefix,
        dst_perc_thr=args.dst_perc_thr,
        edge_len=args.edge_len,
        edge_perc_alr_thr=args.edge_perc_alr_thr,
        max_alr_len_thr=args.max_alr_len_thr,
    )


if __name__ == "__main__":
    raise SystemExit(main())
