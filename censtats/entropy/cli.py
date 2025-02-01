import os
import math
import argparse
import polars as pl
import matplotlib.pyplot as plt
from typing import Generator, TextIO, TYPE_CHECKING, Any
from collections import Counter
# from concurrent.futures import ProcessPoolExecutor

from matplotlib.colors import LinearSegmentedColormap
from intervaltree import Interval, IntervalTree
from loguru import logger

from .constants import DEF_FILTER_RP, DEF_BED9_COLS, DEF_WINDOW_SIZE
from ..common import merge_itvs

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def get_shannon_index_itvs(
    df: pl.DataFrame,
    window_size: int = DEF_WINDOW_SIZE,
    filter_repeats: set[str] | None = None,
) -> Generator[Interval, None, None]:
    """
    Calculate windowed shannon index from repeat content.

    # Args
    * df
            * Dataframe with columns: `["chromStart", "chromEnd", and "name"]`
    * window_size
            * Window size in bases to calculate shannon index over.
            * By default, 5000 bp.
    * filter_repeats
            * Repeats to filter out.
            * By default, `SAR` and `Simple_repeat`

    # Returns
    Generator of `Interval`s with shannon index in `data` attribute.
    """
    if not filter_repeats:
        filter_repeats = DEF_FILTER_RP

    complete_interval = Interval(df["chromStart"].min(), df["chromEnd"].max())
    windows = list(range(complete_interval.begin, complete_interval.end, window_size))
    intervals_windows = IntervalTree()
    for i, start in enumerate(windows):
        try:
            stop = windows[i + 1]
        except IndexError:
            stop = df["chromEnd"].max()
        intervals_windows.addi(start, stop)

    intervals_df = IntervalTree.from_tuples(
        (row["chromStart"], row["chromEnd"], row["name"])
        for row in df.filter(~pl.col("name").is_in(filter_repeats)).iter_rows(
            named=True
        )
    )
    intervals_clipped_windows = []

    for wd in sorted(intervals_windows.iter()):
        overlap = intervals_df.overlap(wd)
        intervals_overlap: IntervalTree = IntervalTree(overlap)
        # Subtract left and right.
        intervals_overlap.chop(complete_interval.begin, wd.begin)
        intervals_overlap.chop(wd.end, complete_interval.end)

        intervals_clipped_windows.append(intervals_overlap)

    for intervals_window in intervals_clipped_windows:
        rp_cnts: Counter[str] = Counter()
        if intervals_windows.is_empty():
            continue
        window_st = intervals_window.begin()
        window_end = intervals_window.end()

        if window_st == 0:
            continue

        for interval_rp in intervals_window:
            rp_cnts[interval_rp.data] += interval_rp.length()
        rp_prop = [cnt / rp_cnts.total() for _, cnt in rp_cnts.items()]
        # https://www.statology.org/shannon-diversity-index/
        num_rp = len(rp_cnts)
        if num_rp <= 1:
            sh_idx = 0.0
        else:
            sh_entropy = -sum(p * math.log(p) for p in rp_prop)
            sh_idx = float(sh_entropy / math.log(num_rp))

        yield Interval(window_st, window_end, round(sh_idx, 3))


def calculate_single_windowed_shannon_index(
    df_group: tuple[tuple[object, ...], pl.DataFrame],
    outdir: str | None,
    window: int,
    ignore_repeats: list[str],
    *,
    omit_plot: bool,
) -> pl.DataFrame | None:
    """
    Calculate windowed shannon index for a chrom `DataFrame` group.

    # Args
    * df_group
            * DataFrame group from `DataFrame.partition()` or `DataFrame.group_by()`
    * outdir
            * Output directory. If `None`, return `DataFrame`
    * window
            * Window size
    * ignore_regions
            * Regions to ignore
    * omit_plot
            * Do not generate plots.

    # Returns
    `DataFrame` if no `outdir`.
    """
    grp, df = df_group
    chrom = grp[0]
    logger.info(f"Calculating Shannon index for {chrom}")
    itvs = list(get_shannon_index_itvs(df, window, filter_repeats=set(ignore_repeats)))
    merged_itvs = merge_itvs(
        itvs,
        dst=1,
        # Merge if index equal.
        fn_cmp=lambda x, y: x.data == y.data,
        fn_merge_itv=lambda x, y: Interval(x.begin, y.end, x.data),
    )
    # Scale colors based on index
    cmap = LinearSegmentedColormap.from_list("", ["red", "orange", "green"])
    df_entropy = pl.DataFrame(
        [
            (
                chrom,
                i.begin,
                i.end,
                "shannon_index",
                i.data,
                "+",
                i.begin,
                i.end,
                # Convert scaled color to rgb
                ",".join(str(round(clr * 255)) for clr in cmap(i.data)[0:-1]),
            )
            for i in merged_itvs
        ],
        schema=DEF_BED9_COLS,
        orient="row",
    )
    if not outdir:
        return df_entropy

    if not omit_plot:
        logger.info(f"Generating plot for {chrom}")
        plt.clf()
        _, ax = plt.subplots(figsize=(20, 10))

        # Draw repeat colors
        for r in df.iter_rows(named=True):
            if r["name"] == "ALR/Alpha" and not r.get("rgb"):
                color = "red"
            else:
                color = r.get("rgb")
            ax.axvspan(
                xmin=r["chromStart"], xmax=r["chromEnd"], facecolor=color, alpha=0.3
            )

        # Plot original values.
        begin, entropy = zip(*[(itv.begin, itv.data) for itv in itvs])
        ax.plot(begin, entropy, color="black")
        ax.fill_between(begin, entropy, color="black")
        ax.margins(x=0, y=0)

        plt.title(f"{chrom} ({window=:,}bp)")
        plt.xlabel("Position")
        plt.ylabel("Shannon index")
        plt.minorticks_on()
        plt.savefig(os.path.join(outdir, f"{chrom}.png"), bbox_inches="tight")
        plt.close()

    df_entropy.write_csv(
        os.path.join(outdir, f"{chrom}.bed"), separator="\t", include_header=False
    )
    return None


def calculate_windowed_shannon_index(
    infile: TextIO,
    outdir: str,
    window_size: int,
    ignore_repeats: list[str],
    # cores: int,
    *,
    omit_plot: bool,
) -> int:
    df_all = pl.read_csv(infile, separator="\t", has_header=False)
    n_cols = df_all.shape[1]
    df_all = df_all.rename(
        dict(zip([f"column_{c}" for c in range(1, n_cols)], DEF_BED9_COLS[:n_cols]))
    )
    os.makedirs(outdir, exist_ok=True)

    df_grps = df_all.partition_by(["chrom"], as_dict=True, maintain_order=True).items()
    for grp in df_grps:
        calculate_single_windowed_shannon_index(
            grp, outdir, window_size, ignore_repeats, omit_plot=omit_plot
        )

    # with ProcessPoolExecutor(cores) as pool:
    # _ = pool.map(
    #     calculate_single_windowed_shannon_index,
    #     *zip(
    #         *[
    #             (grp, outdir, window_size, ignore_regions)
    #             for grp in df_all.partition_by(["chrom"], as_dict=True, maintain_order=True).items()
    #         ]
    #     )
    # )
    return 0


def add_entropy_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "entropy",
        description="Calculate shannon index across a region from RepeatMasker repeats.",
    )
    ap.add_argument(
        "-i",
        "--input",
        required=True,
        help=f"Repeatmasker bedfile. Expects {DEF_BED9_COLS}. The 'name' column should correspond to the repeat name.",
    )
    ap.add_argument(
        "-w",
        "--window",
        type=int,
        default=DEF_WINDOW_SIZE,
        help=f"Window size. Default: {DEF_WINDOW_SIZE}",
    )
    ap.add_argument(
        "-o",
        "--outdir",
        type=str,
        required=True,
        help=(
            "Output dir. Will produce a BED9 file where 'score' corresponds to the Shannon index. "
            "The plot visualize this index across the given repeat region."
        ),
    )
    # ap.add_argument("-c", "--cores", type=int, default=4, help="Number of cores to use.")
    ap.add_argument(
        "--ignore_repeats",
        nargs="*",
        default=[],
        help=f"Repeat types to ignore in calculation. Default: {DEF_FILTER_RP}",
    )
    ap.add_argument("--omit_plot", action="store_true", help="Omit plot.")

    return None
