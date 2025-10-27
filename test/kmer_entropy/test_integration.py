import os
import pytest

from ..helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["input_fa", "output_bed", "expected_bed", "args"],
    [
        (
            f"test/kmer_entropy/input/{chrom}.fa",
            f"test/kmer_entropy/expected/{chrom}.bed",
            f"test/kmer_entropy/expected/{chrom}_expected.bed",
            tuple(["-w", str(100_000), "-k", str(100), "-c", str(4), "--omit_plot"]),
        )
        for chrom in (
            "chr8:40000000-45000000",
            "chr8:40000000-60000000",
        )
    ],
)
def test_check_shannon_entropy(
    input_fa: str, output_bed: str, expected_bed: str, args: tuple[str, ...]
):
    run_integration_test(
        "python",
        "-m",
        "censtats.main",
        "kmer-entropy",
        "-i",
        input_fa,
        "-o",
        os.path.dirname(output_bed),
        *args,
        cmd_output=output_bed,
        expected_output=expected_bed,
    )
