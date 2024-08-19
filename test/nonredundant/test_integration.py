import pytest

from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    [
        "input_left",
        "input_right",
        "expected_output_left",
        "expected_output_right",
        "expected_output_both",
        "expected_output_unaccounted_left",
        "expected_output_unaccounted_right",
        "expected_output_dupe_left",
        "expected_output_dupe_right",
    ],
    [
        (
            "test/nonredundant/input/verkko_all_AS-HOR_lengths.tsv",
            "test/nonredundant/input/hifiasm_all_AS-HOR_lengths.tsv",
            "test/nonredundant/expected/uniq_left.tsv",
            "test/nonredundant/expected/uniq_right.tsv",
            "test/nonredundant/expected/both.tsv",
            "test/nonredundant/expected/unaccounted_left.tsv",
            "test/nonredundant/expected/unaccounted_right.tsv",
            "test/nonredundant/expected/dupes_left.tsv",
            "test/nonredundant/expected/dupes_right.tsv",
        ),
    ],
)
def test_check_cens_status(
    input_left: str,
    input_right: str,
    expected_output_left: str,
    expected_output_right: str,
    expected_output_both: str,
    expected_output_unaccounted_left: str,
    expected_output_unaccounted_right: str,
    expected_output_dupe_left: str,
    expected_output_dupe_right: str,
):
    run_integration_test(
        "python",
        "-m",
        "censtats.main",
        "nonredundant",
        "-l",
        input_left,
        "-r",
        input_right,
        expected_output=[
            ("-ol", expected_output_left),
            ("-or", expected_output_right),
            ("-ul", expected_output_unaccounted_left),
            ("-ur", expected_output_unaccounted_right),
            ("-dl", expected_output_dupe_left),
            ("-dr", expected_output_dupe_right),
            ("-b", expected_output_both),
        ],
    )
