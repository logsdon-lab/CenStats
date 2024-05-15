import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["input_rm_out", "expected_rc_list", "additional_args"],
    [
        (
            "test/status/input/chr21_cens.fa.out",
            "test/status/expected/correct_chr21_cens.tsv",
            (),
        ),
        (
            "test/status/input/chr22_cens.fa.out",
            "test/status/expected/correct_chr22_cens.tsv",
            (),
        ),
        (
            "test/status/input/chr21_chr13_cens_mismap.fa.out",
            "test/status/expected/correct_chr21_chr13_cens_mismap.tsv",
            tuple(["--restrict_13_21"]),
        ),
        (
            "test/status/input/chr9_cens_partials.fa.out",
            "test/status/expected/correct_chr9_cens_partials.tsv",
            (),
        ),
        (
            "test/status/input/chr4_cens_partials.fa.out",
            "test/status/expected/correct_chr4_cens_partials.tsv",
            (),
        ),
        (
            "test/status/input/chr21_cens_false_neg_mismap.fa.out",
            "test/status/expected/correct_chr21_cens_false_neg_mismap.tsv",
            tuple(["--restrict_13_21"]),
        ),
    ],
)
def test_check_cens_status(
    input_rm_out: str, expected_rc_list: str, additional_args: tuple[str]
):
    run_integration_test(
        "python",
        "-m",
        "censtats.main",
        "status",
        "-i",
        input_rm_out,
        "-r",
        "test/status/input/chm13_chm1_cens_v21.trimmed.fa.noheader.out",
        *additional_args,
        expected_output=expected_rc_list,
    )
