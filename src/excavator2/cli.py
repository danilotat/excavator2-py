"""CLI contract for the port. Commands never pretend to run unfinished stages."""

import argparse
from collections.abc import Sequence
from pathlib import Path

from . import __version__


def positive_int(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return number


def mapq_value(value: str) -> int:
    number = int(value)
    if not 0 <= number <= 255:
        raise argparse.ArgumentTypeError("must be between 0 and 255")
    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="excavator2",
        description="Python/C++ port scaffold. Scientific stages are not implemented yet.",
    )
    parser.add_argument("--version", action="version", version=f"excavator2 {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)
    target = commands.add_parser("target", help="Initialize target regions (not implemented)")
    target.add_argument("--config", "--settings", "-s", type=Path, required=True)

    prepare = commands.add_parser("prepare", help="Prepare read counts (not implemented)")
    prepare.add_argument("--mapq", "-q", type=mapq_value, default=20)

    analyze = commands.add_parser("analyze", help="Call CNVs (not implemented)")
    analyze.add_argument("--input", "-i", type=Path, required=True)
    analyze.add_argument(
        "--experiment", "-e", choices=("paired", "pooling", "pooled"), required=True
    )
    analyze.add_argument("--parameters", "-p", type=Path)

    for command in (prepare, analyze):
        command.add_argument("--samples", "-s", type=Path, required=True)
        command.add_argument("--target", "-t", type=Path, required=True)
        command.add_argument("--threads", "-@", type=positive_int, default=1)

    for command in (target, prepare, analyze):
        command.add_argument("--output", "-o", type=Path, required=True)
        command.add_argument("--force", "-f", action="store_true")
        command.add_argument("--verbose", "-v", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    parser.exit(
        status=1,
        message=(
            f"excavator2 {args.command}: this package is a development scaffold; "
            "the scientific stage is not implemented. No output was written. "
            "Use the original scripts for analysis.\n"
        ),
    )
