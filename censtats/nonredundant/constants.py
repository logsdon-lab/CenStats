CHR = {*[f"chr{i}" for i in range(1, 23)], "chrX", "chrY"}
BP_DIFF = 1000
JOIN_COLS = ["sample", "chr", "hap", "ctg", "start", "end", "length"]
JOIN_COLS_RIGHT = [f"{c}_right" for c in JOIN_COLS]
IO_COLS = ["ctg", "start", "end", "length"]
IO_COLS_RIGHT = [f"{c}_right" for c in IO_COLS]
