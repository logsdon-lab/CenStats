import re


RGX_CHR = re.compile(r"chr[0-9XY]{1,2}")
HOR_BP_LEN = 170
DEF_ARR_LEN_THR = 30_000
DEF_BP_JUMP_LEN_THR = 100_000
DEF_EXP_STV_ROW_BED_COLS = ["chr", "start", "stop", "hor", "other", "strand"]
DEF_OUTPUT_BED_COLS = ["chr_name", "start_pos", "stop_pos", "len"]
