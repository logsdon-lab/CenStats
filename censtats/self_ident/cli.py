import os
import pysam
import argparse

from loguru import logger
from typing import TYPE_CHECKING, Any, Generator
from enum import StrEnum
from statistics import mean
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

from .estimate_identity import convertMatrixToBed, createSelfMatrix
from .read_fasta import generateKmersFromFasta


if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


class Dim(StrEnum):
    ONE = "1D"
    TWO = "2D"


def convert_2D_to_1D_ident(
    bed: list[tuple[str, int, int, str, int, int, float]],
    window: int,
    n_bins: int,
    ignore_bands: int,
) -> Generator[tuple[int, int, float], None, None]:
    # Use dictionary to avoid sparse mtx.
    aln_mtx: defaultdict[int, dict[int, float]] = defaultdict(dict)
    for line in bed:
        q_ctg, q_st, q_end, r_ctg, r_st, r_end, perID = line
        # Convert position to indices.
        x = q_st // window
        y = r_st // window
        aln_mtx[x][y] = perID

    st_idxs = list(aln_mtx.keys())
    for st_idx in st_idxs:
        st = st_idx * window + 1
        end = st + window - 1
        band_end_idx = st_idx + n_bins
        # Within the alignment matrix with a n_bins of 5 and ignore_bands of 2:
        # - '*' is the calculated aln band
        # - '+' is self aln.
        # 4 * * *   +
        # 3 * *   +
        # 2 *   +
        # 1   +
        # 0 +
        #   0 1 2 3 4
        mean_ident = mean(
            aln_mtx[x].get(y, 0.0)
            for x in range(st_idx, band_end_idx)
            for y in range(x + ignore_bands, band_end_idx)
        )
        yield st, end, mean_ident


def get_single_self_seq_ident(
    seq_id: str,
    seq: str,
    outdir: str,
    window: int,
    delta: float,
    kmer_size: int,
    ident_thr: float,
    modimizer: int,
    n_bins: int,
    ignore_bands: int,
    dim: Dim,
    round_ndigits: int | None,
) -> None:
    logger.info(f"Generating self sequence identity for {seq_id}.")
    kmers = [kmer_hash for kmer_hash in generateKmersFromFasta(seq, kmer_size)]
    mtx = createSelfMatrix(kmers, window, delta, kmer_size, ident_thr, False, modimizer)
    bed = convertMatrixToBed(mtx, window, ident_thr, seq_id, seq_id, True)
    outfile = os.path.join(outdir, f"{seq_id}.bed")

    if dim == Dim.TWO:
        logger.info(
            f"Writing 2D self sequence identity array for {seq_id} to {outfile}"
        )
        with open(outfile, "wt") as fh:
            for qname, qst, qend, rname, rst, rend, ident in bed:
                ident = round(ident, round_ndigits) if round_ndigits else ident
                fh.write(f"{qname}\t{qst}\t{qend}\t{rname}\t{rst}\t{rend}\t{ident}\n")
    else:
        logger.info(f"Converting 2D self sequence identity matrix to 1D for {seq_id}.")
        logger.info(
            f"Writing 1D self sequence identity array for {seq_id} to {outfile}"
        )
        with open(outfile, "wt") as fh:
            for st, end, ident in convert_2D_to_1D_ident(
                bed, window, n_bins, ignore_bands
            ):
                ident = round(ident, round_ndigits) if round_ndigits else ident
                fh.write(f"{seq_id}\t{st}\t{end}\t{ident}\n")

    return None


def add_self_ident_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "self-ident",
        description="Approximate 1D or 2D average self nucleotide identity via a k-mer-based containment index. Uses ModDotPlot's library.",
    )
    ap.add_argument(
        "-i",
        "--infile",
        required=True,
        type=str,
        help="Input fasta.",
    )
    ap.add_argument(
        "-o",
        "--outdir",
        type=str,
        required=True,
        help="Output directory for self-identity alignment bedfile by contig.",
    )
    ap.add_argument(
        "-x",
        "--dim",
        choices=[Dim.ONE, Dim.TWO],
        help="Dimensionality of self-identity returned.",
        default=Dim.ONE,
        type=Dim,
    )
    ap.add_argument(
        "-p", "--processes", default=4, type=int, help="Number of processes."
    )
    ap.add_argument(
        "-t", "--ident_thr", default=0.86, type=float, help="Identity threshold."
    )
    ap.add_argument("-w", "--window", default=5000, type=int, help="Window size.")
    ap.add_argument("-k", "--kmer_size", default=21, type=int, help="K-mer size.")
    ap.add_argument(
        "-d",
        "--delta",
        default=0.5,
        type=float,
        help="Fraction of neighboring partition to include in identity estimation. Must be between 0 and 1, use > 0.5 is not recommended.",
    )
    ap.add_argument(
        "-m",
        "--modimizer",
        default=1000,
        type=int,
        help="Modimizer sketch size. A lower value will reduce the number of modimizers, but will increase performance. Must be less than --window.",
    )
    # 1D params
    ap.add_argument(
        "-b",
        "--n_bins",
        default=5,
        type=int,
        help="Number of bins to calculate average sequence identity over. Only applicable if mode is 1D.",
    )
    ap.add_argument(
        "--ignore_bands",
        default=2,
        type=int,
        help="Number of bands ignored along self-identity diagonal. Only applicable if mode is 1D.",
    )
    ap.add_argument(
        "--round_ndigits",
        default=None,
        type=int,
        help="Round identity to specified ndigits.",
    )
    return None


def get_self_seq_ident(
    infile: str,
    outdir: str,
    window: int,
    delta: float,
    kmer_size: int,
    ident_thr: float,
    modimizer: int,
    n_bins: int,
    ignore_bands: int,
    processes: int,
    dim: Dim,
    round_ndigits: int | None,
):
    os.makedirs(outdir, exist_ok=True)

    seq = pysam.FastaFile(infile)
    with ProcessPoolExecutor(max_workers=processes) as pool:
        _ = pool.map(
            get_single_self_seq_ident,
            *zip(
                *[
                    (
                        seq_id,
                        seq.fetch(seq_id),
                        outdir,
                        window,
                        delta,
                        kmer_size,
                        ident_thr,
                        modimizer,
                        n_bins,
                        ignore_bands,
                        dim,
                        round_ndigits,
                    )
                    for seq_id in seq.references
                ]
            ),
        )
