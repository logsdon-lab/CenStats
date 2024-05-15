import pytest

from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["input_stv_row_bed", "expected_hor_len_bed"],
    [
        (
            "test/length/input/AS-HOR-vs-chm13_cens_v18.stv_row.all.bed",
            "test/length/expected/AS-HOR_chm13_lengths.bed",
        ),
    ],
)
def test_check_cens_status(
    input_stv_row_bed: str,
    expected_hor_len_bed: str,
):
    run_integration_test(
        "python",
        "-m",
        "censtats.main",
        "length",
        "-i",
        input_stv_row_bed,
        expected_output=expected_hor_len_bed,
    )
