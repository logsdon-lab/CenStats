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
ACROCENTRIC_CHROMOSOMES = {
    "chr13",
    "chr14",
    "chr15",
    "chr21",
    "chr22",
}
RGX_CHR = re.compile(r"(chr[0-9XY]+)")
DST_PERC_THR = 0.3
EDGE_LEN = 100_000
EDGE_PERC_ALR_THR = 0.7
HOR_LEN_THR = 200_000
