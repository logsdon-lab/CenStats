from enum import StrEnum, auto
import polars as pl


class Side(StrEnum):
    Left = auto()
    Right = auto()


BP_DIFF = 1000
IO_COLS = ["ctg", "start", "end", "length"]


def select_exprs(side: Side = Side.Left) -> dict[str, pl.Expr]:
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
