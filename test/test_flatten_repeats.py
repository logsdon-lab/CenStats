import pytest
import polars as pl
from censtats.acrocentrics import flatten_repeats
from censtats.reader import read_repeatmasker_output


@pytest.mark.parametrize(
    ["input_rm", "contig", "df_exp"],
    [
        (
            "test/input/chr22_cens.fa.out",
            "HG00171_chr22_haplotype1-0000002:2413947-12672095",
            pl.DataFrame(
                {"start": 4647309, "end": 8633895, "type": "ALR/Alpha", "dst": 3986586}
            ),
        ),
        (
            "test/input/chr21_cens.fa.out",
            "HG00171_chr21_haplotype1-0000029:4196961-10792258",
            pl.DataFrame(
                {"start": 4560509, "end": 6086655, "type": "ALR/Alpha", "dst": 1526146}
            ),
        ),
    ],
)
def test_flatten_repeats_to_find_largest_alr(input_rm, contig, df_exp):
    df_all = read_repeatmasker_output(input_rm).collect()
    df_ctg = df_all.filter(pl.col("contig") == contig)
    max_alr = (
        flatten_repeats(df_ctg)
        .filter(pl.col("type") == "ALR/Alpha")
        .filter(pl.col("dst") == pl.col("dst").max())
    )
    assert max_alr.equals(df_exp)


@pytest.mark.parametrize(
    ["input_rm", "contig", "df_exp"],
    [
        (
            "test/input/chr21_cens.fa.out",
            "HG00171_chr21_haplotype2-0000154:33297526-38584922",
            pl.DataFrame(
                {
                    "start": [467742, 3990526, 4568967],
                    "end": [2074023, 4224328, 4778825],
                    "type": ["ALR/Alpha"] * 3,
                    "dst": [1606281, 233802, 209858],
                }
            ),
        ),
        (
            "test/input/chr22_cens.fa.out",
            "chm13_chr22:8000001-17400000",
            pl.DataFrame(
                {
                    "start": [486583, 1944302, 4763799],
                    "end": [814885, 2188990, 7786187],
                    "type": ["ALR/Alpha"] * 3,
                    "dst": [328302, 244688, 3022388],
                }
            ),
        ),
    ],
)
def test_flatten_repeats_find_all_major_alrs(input_rm, contig, df_exp):
    df_all = read_repeatmasker_output(input_rm).collect()
    df_ctg = df_all.filter(pl.col("contig") == contig)
    major_alrs = (
        flatten_repeats(df_ctg)
        .filter(pl.col("type") == "ALR/Alpha")
        .filter(pl.col("dst") > 200_000)
    )
    assert major_alrs.equals(df_exp)
