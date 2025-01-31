import argparse
from typing import Any, TYPE_CHECKING

from .status.cli import add_status_cli, check_cens_status
from .length.cli import add_hor_length_cli, calculate_hor_length
from .nonredundant.cli import add_nonredundant_cli, get_nonredundant_cens
from .entropy.cli import add_entropy_cli, calculate_windowed_shannon_index
from .self_ident.cli import add_self_ident_cli, get_self_seq_ident


if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def main() -> int:
    ap = argparse.ArgumentParser(description="Centromere statistics toolkit.")
    sub_ap = ap.add_subparsers(dest="cmd")
    add_status_cli(sub_ap)
    add_hor_length_cli(sub_ap)
    add_nonredundant_cli(sub_ap)
    add_entropy_cli(sub_ap)
    add_self_ident_cli(sub_ap)

    args = ap.parse_args()

    if args.cmd == "status":
        return check_cens_status(
            args.input,
            args.output,
            args.reference,
            reference_prefix=args.reference_prefix,
            dst_perc_thr=args.dst_perc_thr,
            edge_len=args.edge_len,
            edge_perc_alr_thr=args.edge_perc_alr_thr,
            max_alr_len_thr=args.max_alr_len_thr,
            restrict_13_21=args.restrict_13_21,
            restrict_14_22=args.restrict_14_22,
        )
    elif args.cmd == "length":
        return calculate_hor_length(
            infile=args.input_stv,
            rmfile=args.input_rm,
            bp_merge_units=args.bp_merge_units,
            bp_merge_blks=args.bp_merge_blks,
            min_blk_hor_units=args.min_blk_hor_units,
            min_arr_hor_units=args.min_arr_hor_units,
            min_arr_len=args.min_arr_len,
            min_arr_prop=args.min_arr_prop,
            output=args.output,
            output_strand=args.output_strand,
        )
    elif args.cmd == "nonredundant":
        return get_nonredundant_cens(
            args.infile_left,
            args.infile_right,
            args.outfile_left,
            args.outfile_right,
            args.outfile_both,
            args.duplicates_left,
            args.duplicates_right,
            bp_diff=args.diff_bp,
        )
    elif args.cmd == "entropy":
        return calculate_windowed_shannon_index(
            args.input,
            args.outdir,
            args.window,
            args.ignore_repeats,
            # args.cores,
            omit_plot=args.omit_plot,
        )
    elif args.cmd == "self-ident":
        return get_self_seq_ident(
            args.infile,
            args.outdir,
            args.window,
            args.delta,
            args.kmer_size,
            args.ident_thr,
            args.modimizer,
            args.n_bins,
            args.ignore_bands,
            args.processes,
            args.dim,
            args.round_ndigits,
        )
    else:
        raise ValueError(f"Unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
