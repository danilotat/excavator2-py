"""Pack raw preparation oracle exports, including intentional failing cases."""

import argparse
import json
from pathlib import Path

import numpy as np

from excavator2.artifacts import digest, read_export


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exports", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "baseline": json.loads(Path(__file__).with_name("baseline.json").read_text()),
        "files": {},
    }
    for case in sorted(args.exports.glob("*.RData")):
        arrays = {folder.name: read_export(folder) for folder in case.iterdir()}
        destination = args.output / f"{case.stem}.npz"
        np.savez_compressed(destination, **arrays)
        manifest["files"][destination.name] = digest(destination)
    for name in [
        "tools/oracle/characterize_preparation.R",
        "excavator2/lib/F77/F4R.f",
        "excavator2/lib/R/MakeReadCount.R",
        "excavator2/lib/R/LibraryExomeRC.R",
    ]:
        manifest["files"][name] = digest(name)
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
