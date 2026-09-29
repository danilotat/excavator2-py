"""Convert export.R checkpoints into versioned analysis inputs (no R at runtime)."""

import argparse
from pathlib import Path

from excavator2.artifacts import convert_legacy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["prepared-exports", "target-exports", "chromosomes", "centromeres", "output"]:
        parser.add_argument(f"--{name}", type=Path, required=name != "prepared-exports")
    parser.add_argument("--assembly", required=True)
    args = parser.parse_args()
    convert_legacy(
        args.prepared_exports,
        args.target_exports,
        args.chromosomes,
        args.centromeres,
        args.assembly,
        args.output,
    )


if __name__ == "__main__":
    main()
