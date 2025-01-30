from collections import deque
from typing import Callable, Iterable

from intervaltree import Interval


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
