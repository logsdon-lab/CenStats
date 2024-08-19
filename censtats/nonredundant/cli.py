import argparse
from enum import StrEnum, auto
from typing import Any, TextIO, TYPE_CHECKING

import polars as pl
import polars.selectors as cs
from loguru import logger

from .constants import BP_DIFF, JOIN_COLS, JOIN_COLS_RIGHT, IO_COLS


if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


class Side(StrEnum):
    Left = auto()
    Right = auto()


def add_uniq_id(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        id=pl.col("sample")
        + pl.col("chr")
        + pl.col("length").cast(pl.String)
        + pl.col("ctg"),
        id_right=pl.col("sample_right")
        + pl.col("chr_right")
        + pl.col("length_right").cast(pl.String)
        + pl.col("ctg"),
    )


def select_exprs(side: Side) -> dict[str, pl.Expr]:
    suffix = "_right" if side == Side.Right else ""
    return {
        f"ctg{suffix}": pl.col(f"sample{suffix}")
        + "_"
        + pl.when(pl.col(f"rc{suffix}"))
        .then("rc-" + pl.col(f"chr{suffix}"))
        .otherwise(pl.col(f"chr{suffix}"))
        + "_"
        + pl.col(f"ctg{suffix}"),
        f"start{suffix}": pl.col(f"start{suffix}"),
        f"end{suffix}": pl.col(f"end{suffix}"),
        f"length{suffix}": pl.col(f"length{suffix}"),
    }


def restore_name(df: pl.DataFrame, *, side: Side) -> pl.DataFrame:
    """
    Restore original df name.
    """
    cols = select_exprs(side)
    return df.with_columns(**cols).select(cols.keys())


def read_length_tsv(file: str | TextIO) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Read length TSV file
    * Expects columns: `["ctg", "start", "end", "length"]`
    * The `ctg` column must contain info for `["sample", "chr", "ctg"]` and be `'_'` delimited.
    * `ctg` should start with the haplotype information. Either 1 or 2.

    Example:
    * `HG01114_rc-chr1_h2tg000002l#1-130810013:121319346-129631944`
    * `HG01573_rc-chr1_haplotype1-0000024:121168122-126852171`
    """
    df = pl.read_csv(
        file,
        separator="\t",
        has_header=False,
        new_columns=IO_COLS,
    )
    df_formatted = (
        df.with_columns(pl.col("ctg").str.split_exact("_", 2))
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
    return df, df_formatted


def get_nonredundant_cens(
    infile_left: TextIO,
    infile_right: TextIO,
    outfile_left: str,
    outfile_right: str,
    outfile_both: str,
    outfile_unaccounted_left: str,
    outfile_unaccounted_right: str,
    outfile_dupe_left: str,
    outfile_dupe_right: str,
    *,
    bp_diff: int = BP_DIFF,
):
    # Read AS-HOR length dataframe.
    # Calculate cumulative AS-HOR array length per centromere
    # Parse haplotype, chr, sample, and ctg_num_coord.
    df_left_og, df_left_og_fmt = read_length_tsv(infile_left)
    df_right_og, df_right_og_fmt = read_length_tsv(infile_right)

    # Define checks.
    expr_check_left_only = pl.col("length_right").is_null()
    expr_check_right_only = pl.col("length").is_null()
    expr_check_same = pl.col("length_diff").abs() < bp_diff

    # Perform outer join on sample, chr, and hap.
    # Haplotype mappings between assemblies are inconsistent but we can use cumulative HOR array length to infer pairs.
    df_contigs_outer_join = df_left_og_fmt.join(
        df_right_og_fmt, on=["sample", "chr"], how="full"
    ).with_columns(length_diff=(pl.col("length") - pl.col("length_right")).abs())

    # Filter starting dataframes based on checks.
    df_both = df_contigs_outer_join.filter(expr_check_same)

    # Iterate thru the unclear/both mappings and attempt to resolve them.
    unclear_rows_to_remove = []
    left_rows_to_add = []
    right_rows_to_add = []

    # Full join causes duplicates that need to be removed
    rows_both = []
    for _, df_grp in df_both.group_by(["sample", "chr"]):
        df_grp = df_grp.with_row_index()
        used_left_pairs = set()
        used_right_pairs = set()
        # Here we iterate through groups row-by-row and sorts them into each category (left, right, or both).
        while not df_grp.is_empty():
            grp_closest = df_grp.filter(
                pl.col("length_diff") == pl.col("length_diff").min()
            ).row(0, named=True)
            # Generate unique ID for pairs.
            pair_left = (
                grp_closest["chr"] + grp_closest["hap"] + str(grp_closest["length"])
            )
            pair_right = (
                grp_closest["chr_right"]
                + grp_closest["hap_right"]
                + str(grp_closest["length_right"])
            )
            if pair_left not in used_left_pairs and pair_right not in used_right_pairs:
                rows_both.append(grp_closest)
                used_left_pairs.add(pair_left)
                used_right_pairs.add(pair_right)
            elif pair_left not in used_left_pairs and pair_right in used_right_pairs:
                left_rows_to_add.append(grp_closest)
                used_left_pairs.add(pair_left)
            elif pair_left in used_left_pairs and pair_right not in used_right_pairs:
                right_rows_to_add.append(grp_closest)
                used_right_pairs.add(pair_right)
            else:
                left_rows_to_add.append(grp_closest)
                right_rows_to_add.append(grp_closest)

            df_grp = df_grp.filter(pl.col("index") != grp_closest["index"])

    df_both = add_uniq_id(pl.DataFrame(rows_both))
    df_unclear = add_uniq_id(df_contigs_outer_join.filter(~expr_check_same))
    df_left_only = add_uniq_id(df_contigs_outer_join.filter(expr_check_left_only))
    df_right_only = add_uniq_id(df_contigs_outer_join.filter(expr_check_right_only))

    for row in df_unclear.iter_rows(named=True):
        same_id = row["id"] in df_both["id"]
        same_id_right = row["id_right"] in df_both["id_right"]

        # Avoid duplicates in separated rows that appear in both.
        if same_id and not same_id_right:
            right_rows_to_add.append(row)
        elif not same_id and same_id_right:
            left_rows_to_add.append(row)
        else:
            # If neither similar, remove unclear row and sort into categories.
            left_rows_to_add.append(row)
            right_rows_to_add.append(row)

        unclear_rows_to_remove.append(row)

    # Unclear df should be empty.
    assert df_unclear.join(
        pl.DataFrame(unclear_rows_to_remove, schema=df_unclear.schema),
        how="anti",
        on=df_unclear.columns,
    ).is_empty(), "Developer error. Unclear df is not empty."

    # Set columns to None to match the right side.
    df_left_only = (
        pl.concat(
            [
                df_left_only,
                pl.DataFrame(left_rows_to_add, schema=df_left_only.schema),
            ]
        )
        .with_columns(
            **{c: None for c in JOIN_COLS_RIGHT},
        )
        .with_columns(length_diff=pl.col("length") - pl.col("length_right"))
        .unique(subset=JOIN_COLS, keep="first")
        .filter(~pl.col("id").is_in(df_both["id"]))
    )
    df_left = restore_name(df_left_only, side=Side.Left)
    df_right_only = (
        pl.concat(
            [
                df_right_only,
                pl.DataFrame(right_rows_to_add, schema=df_right_only.schema),
            ]
        )
        .with_columns(
            **{c: None for c in JOIN_COLS},
        )
        .with_columns(length_diff=pl.col("length") - pl.col("length_right"))
        .unique(subset=JOIN_COLS_RIGHT, keep="first")
        .filter(~pl.col("id_right").is_in(df_both["id_right"]))
    )
    df_right = restore_name(df_right_only, side=Side.Right)
    df_both = pl.concat(
        [restore_name(df_both, side=Side.Left), restore_name(df_both, side=Side.Right)],
        how="horizontal",
    )
    df_left_unaccounted = pl.concat(
        [df_left, df_both.select(~cs.ends_with("right"))]
    ).join(df_left_og, on="ctg", how="anti")
    df_right_unaccounted = (
        pl.concat([df_right, df_both.select(cs.ends_with("right"))])
        .rename({"ctg_right": "ctg"})
        .join(df_right_og, on="ctg", how="anti")
    )

    df_left_potential_dupes = restore_name(
        df_left_only.filter(pl.col("sample").len().over("sample", "chr", "length") > 1),
        side=Side.Left,
    ).sort(by=["length", "ctg"])
    df_right_potential_dupes = restore_name(
        df_right_only.filter(
            pl.col("sample_right")
            .len()
            .over("sample_right", "chr_right", "length_right")
            > 1
        ),
        side=Side.Right,
    ).sort(by=["length_right", "ctg_right"])

    if not df_left_potential_dupes.is_empty():
        logger.warning(
            f"{df_left_potential_dupes.shape[0]} centromeres potentially duplicated in left file."
        )
    if not df_right_potential_dupes.is_empty():
        logger.warning(
            f"{df_right_potential_dupes.shape[0]} centromeres potentially duplicated in right file."
        )

    if not df_left_unaccounted.is_empty():
        logger.warning(
            f"{df_left_unaccounted.shape[0]} centromeres unaccounted for in left file."
        )
    if not df_right_unaccounted.is_empty():
        logger.warning(
            f"{df_right_unaccounted.shape[0]} centromeres unaccounted for in right file."
        )

    logger.info(f"{df_left.shape[0]} unique centromeres in left file.")
    logger.info(f"{df_right.shape[0]} unique centromeres in right file.")
    logger.info(f"{df_both.shape[0]} centromeres shared by both files.")

    df_left.write_csv(outfile_left, include_header=False, separator="\t")
    df_right.write_csv(outfile_right, include_header=False, separator="\t")
    df_both.write_csv(outfile_both, include_header=False, separator="\t")
    df_left_potential_dupes.write_csv(
        outfile_dupe_left, include_header=False, separator="\t"
    )
    df_right_potential_dupes.write_csv(
        outfile_dupe_right, include_header=False, separator="\t"
    )
    df_left_unaccounted.write_csv(
        outfile_unaccounted_left, include_header=False, separator="\t"
    )
    df_right_unaccounted.write_csv(
        outfile_unaccounted_right, include_header=False, separator="\t"
    )


def add_nonredundant_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "nonredundant",
        description="Get non-redundant list of centromeres based on HOR array length.",
    )
    ap.add_argument(
        "-l",
        "--infile_left",
        help=f"Centromere lengths (left). Expects columns: {IO_COLS}. The '{IO_COLS[0]}' column must be able to be split into {JOIN_COLS[0:4]}. ex. ",
        required=True,
    )
    ap.add_argument(
        "-r",
        "--infile_right",
        help="Centromere lengths (right). Same conditions as --infile_left.",
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
        help="Centromeres shared by both files.",
        type=str,
        default="both.tsv",
    )
    ap.add_argument(
        "-ul",
        "--unaccounted_left",
        help="Unaccounted centromeres from --infile_left.",
        type=str,
        default="unaccounted_left.tsv",
    )
    ap.add_argument(
        "-ur",
        "--unaccounted_right",
        help="Unaccounted centromeres from --infile_right.",
        type=str,
        default="unaccounted_right.tsv",
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
        help="Difference in base pair length between two HOR arrays to be considered different.",
        default=BP_DIFF,
    )
