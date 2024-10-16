from typing import Generator
from collections import Counter

import polars as pl
import numpy as np
import matplotlib.pyplot as plt
from intervaltree import Interval, IntervalTree
from scipy.stats import entropy


def calculate_shannon_index(df: pl.DataFrame, window_size: int = 5000, filter_repeats: set[str] = None) -> Generator[Interval, None, None]:
    """
    Expects st, end, and rp columns in df
    """
    if not filter_repeats:
        filter_repeats = {"HSat1A", "Simple_repeat"}

    complete_interval = Interval(df["st"].min(), df["end"].max())
    windows = list(range(complete_interval.begin, complete_interval.end, window_size))
    intervals_windows = IntervalTree()
    for i, start in enumerate(windows):
        try:
            stop = windows[i+1]
        except IndexError:
            stop=df["end"].max()
        intervals_windows.addi(start, stop)

    intervals_df = IntervalTree.from_tuples(
        (row["st"], row["end"], row["rp"])
        for row in df.filter(~pl.col("rp").is_in(filter_repeats)).iter_rows(named=True)
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
        rp_cnts = Counter()
        intervals_window: IntervalTree
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
            sh_idx = 0
        else:
            sh_idx = entropy(rp_prop) / np.log(num_rp)
        
        yield Interval(window_st, window_end, sh_idx)
    
for chrom, df_chr in df_hg731.sort(by="chrom").group_by(["chrom"], maintain_order=True):
    print(chrom[0])
    wd_begins, wd_ent = zip(*[(i.begin, i.data) for i in calculate_shannon_index(df_chr, 30_000)])
    wd_range = np.arange(len(wd_begins))

    plt.clf()
    ax = plt.subplot()
    ax.plot(wd_begins, wd_ent)
    try:
        poly, (resid, rank, sv, rcond) = np.polynomial.Polynomial.fit(wd_range, wd_ent, deg=2, full=True)
        print(poly)
        # Generate new set of values from polynomial and clip to bounds of entropy.
        coef = poly.convert().coef
        vals = np.clip(
            np.polynomial.polynomial.polyval(wd_range, coef),
            min(wd_ent),
            max(wd_ent),
        )

        ss_total = np.sum((wd_ent - np.mean(wd_ent))**2)
        r_squared = 1 - (resid / ss_total)
        print(r_squared)
        ax.plot(wd_begins, vals)
    except Exception as err:
        print(err)
    for r in df_chr.iter_rows(named=True):
        ax.axvspan(xmin=r["st"], xmax=r["end"], facecolor=r["color"], alpha=0.3)
    plt.show()