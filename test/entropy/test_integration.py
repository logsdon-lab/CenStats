import os
import pytest

from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["input_rm_bed", "output_bed", "expected_bed", "args"],
    [
        (
            "test/entropy/input/HG00096_rm.bed",
            "test/entropy/expected/haplotype1-0000001.bed",
            "test/entropy/expected/haplotype1-0000001_expected.bed",
            tuple(["-w", str(100_000), "--omit_plot"]),
        ),
    ],
)
def test_check_shannon_entropy(
    input_rm_bed: str, output_bed, expected_bed: str, args: tuple[str, ...]
):
    run_integration_test(
        "python",
        "-m",
        "censtats.main",
        "entropy",
        "-i",
        input_rm_bed,
        "-o",
        os.path.dirname(output_bed),
        *args,
        cmd_output=output_bed,
        expected_output=expected_bed,
    )
