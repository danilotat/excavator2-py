"""Pack numeric R exports into small, provenance-labelled FastCall fixtures."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exports", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--wrapper", type=Path, required=True)
    parser.add_argument("--segmentation", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "baseline": json.loads(Path(__file__).with_name("baseline.json").read_text()),
        "R_version": "4.4.1",
        "seed": 20260924,
        "source_sha256": sha256(args.source),
        "wrapper_sha256": sha256(args.wrapper),
        "exporter_sha256": sha256(Path(__file__).with_name("export.R")),
        "packer_sha256": sha256(Path(__file__)),
        "segmentation_sha256": sha256(args.segmentation),
        "conversion": "Binary values reshaped in Fortran order, without decimal conversion.",
        "cases": {},
    }
    for case in sorted(args.exports.glob("*.RData")):
        arrays = {}
        for folder in sorted(case.iterdir()):
            metadata = json.loads((folder / "metadata.json").read_text())
            if metadata["type"] not in ("double", "integer", "logical"):
                raise ValueError(f"Not a numeric fixture: {folder}")
            if any(
                np.fromfile(folder / f"{mask}.bin", dtype="<i4").any()
                for mask in ("missing", "nan")
            ):
                raise ValueError(f"Unexpected missing values in finite fixture: {folder}")
            dtype = "<f8" if metadata["type"] == "double" else "<i4"
            values = np.fromfile(folder / "data.bin", dtype=dtype)
            if metadata["dim"] is not None:
                values = values.reshape(metadata["dim"], order="F")
            arrays[folder.name] = values
        destination = args.output / f"{case.stem}.npz"
        np.savez(destination, **arrays)
        manifest["cases"][case.stem] = {
            "sha256": sha256(destination),
            "arrays": {
                name: {"dtype": str(value.dtype), "shape": list(value.shape)}
                for name, value in arrays.items()
            },
        }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
