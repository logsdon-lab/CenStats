import os
import pytest

from ..helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["input_fasta", "output_bed", "expected_bed", "args"],
    [
        # 1D
        (
            "test/self_ident/input/chm13_chr1.fa",
            "test/self_ident/expected/1D/chr1:121119252-127324151.bed",
            "test/self_ident/expected/1D/chr1:121119252-127324151_expected.bed",
            # Round for testing floats.
            tuple(["--round_ndigits", "3"]),
        ),
        # 2D
        (
            "test/self_ident/input/chm13_chr1.fa",
            "test/self_ident/expected/2D/chr1:121119252-127324151.bed",
            "test/self_ident/expected/2D/chr1:121119252-127324151_expected.bed",
            tuple(["--round_ndigits", "3", "--dim", "2D"]),
        ),
    ],
)
def test_check_shannon_entropy(
    input_fasta: str, output_bed: str, expected_bed: str, args: tuple[str, ...]
):
    run_integration_test(
        "python",
        "-m",
        "censtats.main",
        "self-ident",
        "-i",
        input_fasta,
        "-o",
        os.path.dirname(output_bed),
        *args,
        cmd_output=output_bed,
        expected_output=expected_bed,
    )
