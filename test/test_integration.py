import pytest
import subprocess


@pytest.mark.parametrize(
    ["input_rm_out", "expected_rc_list", "additional_args"],
    [
        ("test/input/chr21_cens.fa.out", "test/expected/correct_chr21_cens.tsv", ()),
        # Partially incorrect due to only using q-arm. Actuall chr10.
        # HG00171_chr22_haplotype2-0000170:38603708-44922516.
        ("test/input/chr22_cens.fa.out", "test/expected/correct_chr22_cens.tsv", ()),
        # Partially incorrect. Using only q-arm means repeats to compare. edit distance between fwd and rev close.
        # TODO: fix.
        (
            "test/input/chr21_chr13_cens_mismap.fa.out",
            "test/expected/correct_chr21_chr13_cens_mismap.tsv",
            tuple(["--restrict_13_21"]),
        ),
        (
            "test/input/chr9_cens_partials.fa.out",
            "test/expected/correct_chr9_cens_partials.tsv",
            (),
        ),
        (
            "test/input/chr4_cens_partials.fa.out",
            "test/expected/correct_chr4_cens_partials.tsv",
            (),
        ),
        (
            "test/input/chr21_cens_false_neg_mismap.fa.out",
            "test/expected/correct_chr21_cens_false_neg_mismap.tsv",
            tuple(["--restrict_13_21"]),
        ),
    ],
)
def test_check_cens_status(
    input_rm_out: str, expected_rc_list: str, additional_args: tuple[str]
):
    process = subprocess.run(
        [
            "python",
            "-m",
            "censtats.main",
            "status",
            "-i",
            input_rm_out,
            "-r",
            "test/input/chm13_chm1_cens_v21.trimmed.fa.noheader.out",
            *additional_args,
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
