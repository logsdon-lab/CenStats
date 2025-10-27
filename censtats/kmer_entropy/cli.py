import os
import math
import argparse
import pyfaidx
import polars as pl
import matplotlib.pyplot as plt

from loguru import logger
from collections import Counter
from intervaltree import Interval
from typing import Generator, TYPE_CHECKING, Any
from matplotlib.colors import LinearSegmentedColormap
from concurrent.futures import ProcessPoolExecutor, Future, as_completed

from .constants import DEF_WINDOW_SIZE, DEF_KMER_SIZE, DEF_BED9_COLS
from ..common import generateKmersFromFasta

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def calculate_kmer_shannon_index_itv(qitv: Interval, kmer_size: int, fa: str) -> Interval:
    chrom, start, stop = qitv.data, qitv.begin, qitv.end
    fa_fh = pyfaidx.Fasta(fa)
    seq = fa_fh.faidx.fetch(name=chrom, start=start, end=stop)

    logger.info(f"On {chrom}:{start}-{stop}...")

    window_len = stop - start
    total_kmers = window_len - kmer_size + 1
    kmers = Counter(generateKmersFromFasta(str(seq), kmer_size))

    rp_prop = [cnt / total_kmers for _, cnt in kmers.items()]
    # https://www.statology.org/shannon-diversity-index/
    num_rp = len(kmers)
    if num_rp <= 1:
        sh_idx = 0.0
    else:
        sh_entropy = -sum(p * math.log(p) for p in rp_prop)
        sh_idx = float(sh_entropy / math.log(num_rp))

    return Interval(start, stop, round(sh_idx, 3))


def process_kmer_shannon_index_all_itvs(
    infile: str,
    chrom: str,
    chrom_len: int,
    cores: int,
    window_size: int = DEF_WINDOW_SIZE,
    kmer_size: int = DEF_KMER_SIZE,
) -> Generator[Interval, None, None]:
    """
    Calculate windowed shannon index from repeat content.

    # Args
    * infile
            * Input fasta file.
    * chrom
            * Chromosome name.
    * chrom_len
            * Chromosome length.
    * cores
            * Number of cores.
    * window_size
            * Window size in bases to calculate shannon index over.
            * By default, 5000 bp.
    * kmer_size
            * Kmer size.

    # Returns
    Generator of `Interval`s with shannon index in `data` attribute.
    """
    windows = list(range(1, chrom_len, window_size))

    with ProcessPoolExecutor(max_workers=cores) as pool:
        futures = []
        for i, start in enumerate(windows):
            try:
                stop = windows[i + 1]
            except IndexError:
                stop = chrom_len

            future: Future[Interval] = pool.submit(
                calculate_kmer_shannon_index_itv,
                Interval(start, stop, chrom),
                kmer_size,
                infile,
            )
            futures.append(future)

    for future in as_completed(futures):
        if future.cancelled():
            logger.error(
                f"Failed to calculate shannon index for query interval: {future.exception()}"
            )
            continue
        yield future.result()


def calculate_plot_windowed_kmer_shannon_index(
    infile: str,
    chrom: str,
    chrom_len: int,
    kmer_size: int,
    window_size: int,
    outdir: str | None,
    cores: int,
    *,
    omit_plot: bool,
) -> pl.DataFrame | None:
    """
    Calculate windowed shannon index based on kmer set for a chrom.

    # Args
    * infile
            * Input fasta file.
    * chrom
            * Chromosome name.
    * chrom_len
            * Chromosome length.
    * kmer_size
            * Kmer size.
    * window_size
            * Window size
    * cores
            * Number of cores.
    * outdir
            * Output directory. If `None`, return `DataFrame`
    * omit_plot
            * Do not generate plots.

    # Returns
    `DataFrame` if no `outdir`.
    """

    logger.info(f"Calculating Shannon index for {chrom}")
    itvs = list(
        process_kmer_shannon_index_all_itvs(
            infile=infile,
            chrom=chrom,
            chrom_len=chrom_len,
            window_size=window_size,
            kmer_size=kmer_size,
            cores=cores,
        )
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
            for i in itvs
        ],
        schema=DEF_BED9_COLS,
        orient="row",
    ).sort(by="chromStart")

    if not outdir:
        return df_entropy

    if not omit_plot:
        logger.info(f"Generating plot for {chrom}")
        plt.clf()
        _, ax = plt.subplots(figsize=(10, 2.5))

        # Plot original values.
        if not df_entropy.is_empty():
            begin, entropy = zip(
                *[
                    (itv["chromStart"], itv["score"])
                    for itv in df_entropy.iter_rows(named=True)
                ]
            )
            ax.plot(begin, entropy, color="black")
            ax.fill_between(begin, entropy, color="black")

        ax.margins(x=0, y=0)

        plt.title(f"{chrom} ({window_size=:,} bp, {kmer_size=:,} bp)")
        plt.xlabel("Position")
        plt.ylabel("Shannon index")
        plt.minorticks_on()
        plt.savefig(os.path.join(outdir, f"{chrom}.png"), bbox_inches="tight")
        plt.close()

    df_entropy.write_csv(
        os.path.join(outdir, f"{chrom}.bed"), separator="\t", include_header=False
    )
    return None


def calculate_kmer_windowed_shannon_index(
    infile: str,
    outdir: str,
    window_size: int,
    kmer_size: str,
    cores: int,
    *,
    omit_plot: bool,
) -> int:
    os.makedirs(outdir, exist_ok=True)

    fa = pyfaidx.Fasta(infile)
    fai = fa.faidx

    for chrom, rec in fai.index.items():
        rec: pyfaidx.IndexRecord
        chrom_len = rec.rlen
        calculate_plot_windowed_kmer_shannon_index(
            infile=infile,
            chrom=chrom,
            chrom_len=chrom_len,
            kmer_size=kmer_size,
            window_size=window_size,
            outdir=outdir,
            omit_plot=omit_plot,
            cores=cores,
        )

    return 0


def add_kmer_entropy_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "kmer-entropy",
        description="Calculate shannon index across a region from kmers.",
    )
    ap.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Fasta file of region to evaluate.",
    )
    ap.add_argument(
        "-w",
        "--window",
        type=int,
        default=DEF_WINDOW_SIZE,
        help=f"Window size. Default: {DEF_WINDOW_SIZE}",
    )
    ap.add_argument(
        "-k",
        "--kmer-size",
        type=int,
        default=DEF_KMER_SIZE,
        help=f"Kmer size. Default: {DEF_KMER_SIZE}",
    )
    ap.add_argument(
        "-o",
        "--outdir",
        type=str,
        required=True,
        help=(
            "Output dir. Will produce a BED9 file where 'score' corresponds to the Shannon index. "
            "The plot visualizes this value across the given repeat region."
        ),
    )
    ap.add_argument(
        "-c", "--cores", type=int, default=4, help="Number of cores to use."
    )
    ap.add_argument("--omit_plot", action="store_true", help="Omit plot.")

    return None
