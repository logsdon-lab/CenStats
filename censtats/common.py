from collections import deque
from typing import Callable, Iterable

from intervaltree import Interval

# ModDotPlot
# https://github.com/marbl/ModDotPlot/commit/0f593a7b7b317cdfc00ef350491e17239eda594f
from typing import Generator
import pyfaidx
import mmh3

tab_b = bytes.maketrans(b"ACTG", b"TGAC")


def generateKmersFromFasta(seq: str, k: int) -> Generator[int, None, None]:
    n = len(seq)
    for i in range(n - k + 1):
        # Remove case sensitivity
        kmer = seq[i : i + k].upper()
        fh = mmh3.hash(kmer)

        # Calculate reverse complement hash directly without the need for translation
        rc = mmh3.hash(kmer[::-1].translate(tab_b))

        yield fh if fh < rc else rc


def readKmersFromFile(
    filename: str, ksize: int
) -> Generator[tuple[str, list[int]], None, None]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    seq = pyfaidx.Fasta(filename)

    for seq_rec in seq:
        kmers_for_seq = [
            kmer_hash for kmer_hash in generateKmersFromFasta(str(seq_rec), ksize)
        ]
        yield seq_rec.name, kmers_for_seq


def fn_cmp_def(_itv_1: Interval, _itv_2: Interval) -> bool:
    return True


def fn_merge_itv_def(itv_1: Interval, itv_2: Interval) -> Interval:
    return Interval(begin=itv_1.begin, end=itv_2.end, data=None)


def merge_itvs(
    itvs: Iterable[Interval],
    dst: int = 1,
    fn_cmp: Callable[[Interval, Interval], bool] | None = None,
    fn_merge_itv: Callable[[Interval, Interval], Interval] | None = None,
) -> list[Interval]:
    if not fn_cmp:
        fn_cmp = fn_cmp_def
    if not fn_merge_itv:
        fn_merge_itv = fn_merge_itv_def

    final_itvs = []
    sorted_itvs = deque(sorted(itvs))
    while sorted_itvs:
        try:
            itv_1 = sorted_itvs.popleft()
        except IndexError:
            break
        try:
            itv_2 = sorted_itvs.popleft()
        except IndexError:
            final_itvs.append(itv_1)
            break
        dst_between = itv_2.begin - itv_1.end
        passes_cmp = fn_cmp(itv_1, itv_2)
        if dst_between <= dst and passes_cmp:
            sorted_itvs.appendleft(fn_merge_itv(itv_1, itv_2))
        else:
            final_itvs.append(itv_1)
            sorted_itvs.appendleft(itv_2)

    return final_itvs
