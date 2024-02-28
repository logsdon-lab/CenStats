import polars as pl


def jaccard_index(a: set[str], b: set[str]) -> float:
    """
    Jaccard similarity index.
    * https://www.statisticshowto.com/jaccard-index/
    """
    return (len(a.intersection(b)) / len(a.union(b))) * 100.0


def get_contig_similarity_by_jaccard_index(
    contigs: list[str], ref_contigs: list[str], jaccard_index: list[float]
) -> pl.DataFrame:
    return (
        pl.LazyFrame(
            {"contig": contigs, "ref": ref_contigs, "similarity": jaccard_index}
        )
        .with_columns(
            pl.col("similarity").max().over("contig").alias("highest_similarity")
        )
        .filter(pl.col("similarity") == pl.col("highest_similarity"))
        .select(["contig", "ref", "similarity"])
        .collect()
    )
