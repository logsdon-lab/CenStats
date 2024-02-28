import polars as pl

from .constants import RM_COLS


def read_repeatmasker_output(input_path: str) -> pl.LazyFrame:
    return (
        pl.scan_csv(input_path, separator="\t", has_header=False, new_columns=RM_COLS)
        .with_columns(
            pl.col("type")
            .str.replace_many(
                [
                    "/ERVK",
                    "/ERVL",
                    "/ERV1",
                    "/CR1",
                    "/L1",
                    "/L2",
                    "/RTE-X",
                    "/RTE-BovB",
                    "/Gypsy",
                    "-MaLR",
                    "/Alu",
                    "/Deu",
                    "/MIR",
                    "?",
                    "/hAT",
                    "/hAT-Blackjack",
                    "/hAT-Charlie",
                    "/MULE-MuDR",
                    "/PiggyBac",
                    "/TcMar-Mariner" "/TcMar",
                    "/TcMar?",
                    "/hAT-Tip100" "/TcMar-Tigger" "/Dong-R4" "/tRNA",
                ],
                "",
            )
            .str.replace_many(
                [
                    "DNA-Tc2",
                    "DNA?",
                    "DNA-Blackjack",
                    "DNA-Charlie",
                    "DNA-Tigger",
                    "DNA-Tip100",
                ],
                "DNA",
            )
            .str.replace("GSATX", "GSAT")
            .str.replace("LTR\\S", "LTR")
            .str.replace("SAR", "HSat1A")
            .str.replace("HSAT", "HSat1B")
            .str.replace_many(["HSATII", "(CATTC)n", "(GAATG)n"], "HSat2")
        )
        .with_columns(dst=pl.col("end") - pl.col("start"))
        .drop("div", "deldiv", "insdiv", "x", "y", "z", "left", "right", "idx")
    )
