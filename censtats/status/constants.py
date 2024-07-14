import re


RM_COLS = [
    "idx",
    "div",
    "deldiv",
    "insdiv",
    "contig",
    "start",
    "end",
    "left",
    "C",
    "type",
    "rClass",
    "right",
    "x",
    "y",
    "z",
]
CHROMOSOMES_13_21 = {
    "chr13",
    "chr21",
}
CHROMOSOMES_14_22 = {
    "chr14",
    "chr22",
}

RGX_CHR = re.compile(r"(chr[0-9XY]+)")
DST_PERC_THR = 0.3
EDGE_LEN = 100_000
EDGE_PERC_ALR_THR = 0.7
MAX_ALR_LEN_THR = 200_000
REPEAT_SPLIT_LEN = 2_000
