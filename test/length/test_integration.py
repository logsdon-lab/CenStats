import pytest

from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    [
        "input_stv_row_bed",
        "input_rm_bed",
        "expected_arr_len_bed",
        "expected_arr_len_strand_bed",
        "options",
    ],
    [
        (
            "test/length/input/AS-HOR-vs-chm13_cens_v18.stv_row.all.bed",
            None,
            "test/length/expected/AS-HOR-vs-chm13_cens_v18.stv_row.all.bed",
            "test/length/expected/AS-HOR-vs-chm13_cens_v18.stv_row.all.strand.bed",
            tuple(),
        ),
        # chr10 small. Small HORs.
        (
            "test/length/input/HG03520_chr10_haplotype2-0000161:38555523-44048872.bed",
            None,
            "test/length/expected/HG03520_chr10_haplotype2-0000161:38555523-44048872.bed",
            None,
            tuple(),
        ),
        # chr4 merge across other repeats.
        (
            "test/length/input/NA19036_chr4_haplotype2-0000075:49294287-52381925.bed",
            None,
            "test/length/expected/NA19036_chr4_haplotype2-0000075:49294287-52381925_incorrect.bed",
            None,
            # Large bp_merge_blks value
            tuple(["-mb", "100000"]),
        ),
        # Use repeatmasker tracks to avoid merging across other repeats.
        (
            "test/length/input/NA19036_chr4_haplotype2-0000075:49294287-52381925.bed",
            "test/length/input/NA19036_chr4_haplotype2-0000075:49294287-52381925.out",
            "test/length/expected/NA19036_chr4_haplotype2-0000075:49294287-52381925.bed",
            "test/length/expected/NA19036_chr4_haplotype2-0000075:49294287-52381925_strand.bed",
            # Large bp_merge_blks value
            tuple(["-mb", "100000"]),
        ),
        # chr2 LINE elements.
        (
            "test/length/input/NA19331_chr2_haplotype2-0000182:91704425-96165312.bed",
            "test/length/input/NA19331_chr2_haplotype2-0000182:91704425-96165312.out",
            "test/length/expected/NA19331_chr2_haplotype2-0000182:91704425-96165312.bed",
            "test/length/expected/NA19331_chr2_haplotype2-0000182:91704425-96165312_strand.bed",
            tuple(),
        ),
    ],
)
def test_check_arr_len(
    input_stv_row_bed: str,
    input_rm_bed: str | None,
    expected_arr_len_bed: str,
    expected_arr_len_strand_bed: str | None,
    options: tuple[str, ...],
):
    outputs = [("-o", expected_arr_len_bed)]
    additional_args = [
        *options,
    ]
    if expected_arr_len_strand_bed:
        outputs.append(("-s", expected_arr_len_strand_bed))

    if input_rm_bed:
        additional_args.extend(("-r", input_rm_bed))

    run_integration_test(
        "python",
        "-m",
        "censtats.main",
        "length",
        "-i",
        input_stv_row_bed,
        *additional_args,
        expected_output=outputs,
    )
