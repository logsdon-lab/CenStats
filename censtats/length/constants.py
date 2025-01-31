DEF_MIN_BLK_HOR_UNITS = 2
DEF_MIN_ARR_HOR_UNITS = 10
DEF_BP_MERGE_BLKS = 7000
# 1.5 monomers.
DEF_BP_MERGE_UNITS = 256
DEF_MERGE_RBLACKLIST = {"BSR", "GSAT", "SAR", "HSATII", "HSATI", "(CATTC)n", "(GAATG)n"}

DEF_INPUT_BED_COLS = ["chrom", "chrom_st", "chrom_end", "name", "score", "strand"]
DEF_INPUT_RM_COLS = [
    "contig",
    "start",
    "end",
    "type",
    "rClass",
]
DEF_INPUT_RM_COL_IDX = [4, 5, 6, 9, 10]

DEF_OUTPUT_BED_COLS = ["chrom", "chrom_st", "chrom_end", "name", "score", "prop"]
DEF_OUTPUT_BED_COLS_STRAND = [
    "chrom",
    "chrom_st",
    "chrom_end",
    "name",
    "score",
    "prop",
    "strand",
]
