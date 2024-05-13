import pytest
import subprocess


@pytest.mark.parametrize(
    "input_rm_out,expected_rc_list",
    [
        (
            "test/input/chr21_cens.fa.out",
            "test/expected/correct_chr21_cens.tsv",
        ),
        (
            "test/input/chr22_cens.fa.out",
            "test/expected/correct_chr22_cens.tsv",
        ),
        (
            "test/input/chr21_cens_partials.fa.out",
            "test/expected/correct_chr21_cens_partials.tsv",
        ),
        (
            "test/input/chr21_chr13_cens_mismap.fa.out",
            "test/expected/correct_chr21_chr13_cens_mismap.tsv",
        ),
        (
            "test/input/chr9_cens_partials.fa.out",
            "test/expected/correct_chr9_cens_partials.tsv",
        ),
    ],
)
def test_check_cens_status(input_rm_out: str, expected_rc_list: str):
    process = subprocess.run(
        [
            "censtats",
            "-i",
            input_rm_out,
            "-r",
            "test/input/chm13_chm1_cens_v21.trimmed.fa.noheader.out",
        ],
        capture_output=True,
        check=True,
    )
    res = sorted(
        line.split("\t") for line in process.stdout.decode().split("\n") if line
    )
    with open(expected_rc_list, "rt") as exp_res_fh:
        exp_res = sorted(
            line.strip().split("\t") for line in exp_res_fh.readlines() if line
        )
        assert res == exp_res
