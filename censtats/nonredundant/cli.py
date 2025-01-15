import sys
import argparse
from typing import Any, TextIO, TYPE_CHECKING

import polars as pl
from loguru import logger

from .constants import BP_DIFF, IO_COLS, Side, select_exprs


if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
    RowVals = tuple[Any, ...]
else:
    SubArgumentParser = Any


def read_as_hor_length_tsv(file: str | TextIO) -> pl.DataFrame:
    """
    Read AS-HOR array length TSV file. Groups by ctg and gets aggregated length before formatting.
    * Expects columns: `["ctg", "start", "end", "length"]`
    * The `ctg` column must contain info for `["sample", "chr", "ctg"]` and be `'_'` delimited.
    * `ctg` should start with the haplotype information. Either 1 or 2.

    Example:
    * `HG01114_rc-chr1_h2tg000002l#1-130810013:121319346-129631944`
    * `HG01573_rc-chr1_haplotype1-0000024:121168122-126852171`
    """
    return (
        pl.read_csv(
            file,
            separator="\t",
            has_header=False,
            new_columns=IO_COLS,
        )
        .group_by("ctg")
        .agg(
            pl.col("start").min(), pl.col("end").max(), pl.sum("length").alias("length")
        )
        .with_columns(og_ctg=pl.col("ctg"), ctg=pl.col("ctg").str.split_exact("_", 2))
        .unnest("ctg")
        .rename({"field_0": "sample", "field_1": "chr", "field_2": "ctg"})
        .with_columns(
            # TODO: Handle mat or pat.
            hap=pl.col("ctg")
            .str.extract(r"(^h.*?\d)", 1)
            .replace({"haplotype1": "h1", "haplotype2": "h2"})
            .fill_null("unassigned"),
            rc=pl.col("chr").str.contains("rc-"),
            chr=pl.col("chr").str.replace("rc-", ""),
        )
    )


def get_nonredundant_cens(
    infile_left: TextIO,
    infile_right: TextIO,
    outfile_left: str,
    outfile_right: str,
    outfile_both: str,
    outfile_dupe_left: str,
    outfile_dupe_right: str,
    *,
    bp_diff: int = BP_DIFF,
):
    # Read AS-HOR length dataframe.
    # Calculate cumulative AS-HOR array length per centromere
    # Parse haplotype, chr, sample, and ctg_num_coord.
    df_left_og_fmt = read_as_hor_length_tsv(infile_left)
    df_right_og_fmt = read_as_hor_length_tsv(infile_right)

    df_merged_og_fmt = pl.concat(
        [
            df_left_og_fmt.with_columns(rtype=pl.lit("left")),
            df_right_og_fmt.with_columns(rtype=pl.lit("right")),
        ]
    )
    cens: dict[str, set[RowVals]] = {"left": set(), "right": set()}
    dupe_cens: dict[str, set[RowVals]] = {"left": set(), "right": set()}
    shared_cens: set[tuple[RowVals, RowVals]] = set()
    covered_cens: set[str] = set()

    # Haplotype information is unreliable so have to do row-by-row comparison.
    # Full outer join makes things very difficult.
    for _, df_grp in df_merged_og_fmt.group_by(["sample", "chr"]):
        df_grp = df_grp.filter(~pl.col("og_ctg").is_in(covered_cens))
        for row_1 in df_grp.iter_rows(named=True):
            row_1_vals = tuple(row_1.values())
            for row_2 in df_grp.iter_rows(named=True):
                row_2_vals = tuple(row_2.values())
                if (
                    row_1["og_ctg"] == row_2["og_ctg"]
                    or row_2["og_ctg"] in covered_cens
                ):
                    continue

                abs_bp_diff = abs(row_1["length"] - row_2["length"])
                # Get potential duplicates.
                # Very rare for AS-HOR array length to be exactly identical.
                if abs_bp_diff == 0 and row_1["rtype"] == row_2["rtype"]:
                    dupe_cens[row_1["rtype"]].add(row_1_vals)
                    dupe_cens[row_1["rtype"]].add(row_2_vals)

                if abs_bp_diff > bp_diff:
                    continue

                if row_1["og_ctg"] in covered_cens:
                    continue

                covered_cens.add(row_1["og_ctg"])
                covered_cens.add(row_2["og_ctg"])
                shared_cens.add((row_1_vals, row_2_vals))

            if row_1["og_ctg"] in covered_cens:
                continue

            cens[row_1["rtype"]].add(row_1_vals)
            covered_cens.add(row_1["og_ctg"])

    del covered_cens

    new_schema = df_left_og_fmt.schema
    new_schema["rtype"] = pl.String
    df_left = pl.DataFrame(list(cens["left"]), schema=new_schema, orient="row").select(
        **select_exprs()
    )
    df_right = pl.DataFrame(
        list(cens["right"]), schema=new_schema, orient="row"
    ).select(**select_exprs())

    shared_left, shared_right = zip(*list(shared_cens))
    df_shared = pl.concat(
        [
            pl.DataFrame(shared_left, schema=new_schema, orient="row").select(
                **select_exprs()
            ),
            pl.DataFrame(
                shared_right,
                schema=[f"{s}_right" for s in new_schema.keys()],
                orient="row",
            ).select(**select_exprs(side=Side.Right)),
        ],
        how="horizontal",
    )
    df_left_potential_dupes = (
        pl.DataFrame(list(dupe_cens["left"]), schema=new_schema, orient="row")
        .select(**select_exprs())
        .sort(by=["ctg", "length"])
    )
    df_right_potential_dupes = (
        pl.DataFrame(list(dupe_cens["right"]), schema=new_schema, orient="row")
        .select(**select_exprs())
        .sort(by=["ctg", "length"])
    )

    if not df_left_potential_dupes.is_empty():
        logger.warning(
            f"{df_left_potential_dupes.shape[0]} centromeres potentially duplicated in left file."
        )
    if not df_right_potential_dupes.is_empty():
        logger.warning(
            f"{df_right_potential_dupes.shape[0]} centromeres potentially duplicated in right file."
        )

    logger.info(f"{df_left.shape[0]} unique centromeres in left file.")
    logger.info(f"{df_right.shape[0]} unique centromeres in right file.")
    logger.info(f"{df_shared.shape[0]} centromeres shared by both files.")

    df_left.write_csv(outfile_left, include_header=False, separator="\t")
    df_right.write_csv(outfile_right, include_header=False, separator="\t")
    df_shared.write_csv(outfile_both, include_header=False, separator="\t")
    df_left_potential_dupes.write_csv(
        outfile_dupe_left, include_header=False, separator="\t"
    )
    df_right_potential_dupes.write_csv(
        outfile_dupe_right, include_header=False, separator="\t"
    )


def add_nonredundant_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "nonredundant",
        description="Get non-redundant list of centromeres based on HOR array length from two AS-HOR array length lists.",
    )
    ap.add_argument(
        "-l",
        "--infile_left",
        type=argparse.FileType("rb"),
        help=" ".join(
            [
                f"Centromere lengths. Expects columns: {IO_COLS}.",
                f"The '{IO_COLS[0]}' column must be able to be split by '_' into ['sample', 'chr', 'hap-ctg'].",
                "ex. 'HG01573_rc-chr1_haplotype1-0000024:121168122-126852171'",
            ]
        ),
        required=True,
    )
    ap.add_argument(
        "-r",
        "--infile_right",
        type=argparse.FileType("rb"),
        default=sys.stdin,
        help="Centromere lengths. Same format expected as --infile_left.",
        required=True,
    )
    ap.add_argument(
        "-ol",
        "--outfile_left",
        help="Unique centromeres from --infile_left.",
        type=str,
        default="uniq_left.tsv",
    )
    ap.add_argument(
        "-or",
        "--outfile_right",
        help="Unique centromeres from --infile_right.",
        type=str,
        default="uniq_right.tsv",
    )
    ap.add_argument(
        "-b",
        "--outfile_both",
        help="Centromeres shared by both files. Centromeres from --infile_left in first four columns and centromeres from --infile_right in last four columns.",
        type=str,
        default="both.tsv",
    )
    ap.add_argument(
        "-dl",
        "--duplicates_left",
        help="Potentially duplicated centromeres from --infile_left.",
        type=str,
        default="dupes_left.tsv",
    )
    ap.add_argument(
        "-dr",
        "--duplicates_right",
        help="Potentially duplicated centromeres from --infile_right.",
        type=str,
        default="dupes_right.tsv",
    )
    ap.add_argument(
        "-d",
        "--diff_bp",
        type=int,
        help=f"Difference in base pair length between two HOR arrays to be considered different. Defaults to {BP_DIFF} bp",
        default=BP_DIFF,
    )
