import argparse
from typing import Any, TYPE_CHECKING

from .status.cli import add_status_cli, check_cens_status

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def main() -> int:
    ap = argparse.ArgumentParser(description="Centromere statistics tool kit.")
    sub_ap = ap.add_subparsers(dest="cmd")
    add_status_cli(sub_ap)

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
        raise NotImplementedError("Length command not implemented.")
    else:
        raise ValueError(f"Unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
