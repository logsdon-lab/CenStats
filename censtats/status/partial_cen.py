import polars as pl

from .constants import EDGE_LEN, EDGE_PERC_ALR_THR, MAX_ALR_LEN_THR


def is_partial_centromere(
    df: pl.DataFrame,
    *,
    edge_len: int = EDGE_LEN,
    edge_perc_alr_thr: float = EDGE_PERC_ALR_THR,
    max_alr_len_thr: int = MAX_ALR_LEN_THR,
) -> bool:
    """
    Check if centromere is partially constructed based on ALR percentage at either ends of the contig.

    ALR/Alpha repeat content is used as a heuristic as we expect primarily monomeric repeats to be at the edges of the HOR array.

    ### Args
    `df`
        RepeatMasker output with a `dst` column for each repeat row.
    `edge_len`
        Edge len to check. Defaults to 100 kbp.
    `edge_perc_alr_thr`
        ALR percentage threshold needed to be considered incomplete. Defaults to 70%.
    `max_alr_len_thr`
        Length threshold of largest ALR needed to be considered complete. Defaults to 200 kbp.

    ### Returns
    Whether the centromere is partially constructed.
    """
    # Check if partial centromere based on ALR perc on ends.
    # Check N kbp from start and end of contig.
    ledge = (
        df.filter(pl.col("start") < edge_len)
        # Ensure that don't get more than edge length.
        .with_columns(
            end=pl.when(pl.col("end") > edge_len)
            .then(pl.lit(edge_len))
            .otherwise(pl.col("end"))
        )
        .with_columns(dst=pl.col("end") - pl.col("start"))
    )
    redge = df.filter(pl.col("start") > df[-1]["end"] - edge_len)
    try:
        ledge_perc_alr = (
            ledge.group_by("type")
            .agg(pl.col("dst").sum() / ledge["dst"].sum())
            .filter(pl.col("type") == "ALR/Alpha")
            .row(0)[1]
        )
    except pl.exceptions.OutOfBoundsError:
        ledge_perc_alr = 0.0 if not ledge.is_empty() else 100.0
    try:
        redge_perc_alr = (
            redge.group_by("type")
            .agg(pl.col("dst").sum() / redge["dst"].sum())
            .filter(pl.col("type") == "ALR/Alpha")
            .row(0)[1]
        )
    except pl.exceptions.OutOfBoundsError:
        redge_perc_alr = 0.0 if not redge.is_empty() else 100.0

    # Check if edges have ALR.
    are_edges_alr = (
        ledge_perc_alr > edge_perc_alr_thr or redge_perc_alr > edge_perc_alr_thr
    )
    # If they don't, check that the contig has at least one ALR that meets the minimum threshold len.
    if not are_edges_alr:
        max_alr_len = df.filter(pl.col("type") == "ALR/Alpha").get_column("dst").max()
        return max_alr_len < max_alr_len_thr

    return are_edges_alr
