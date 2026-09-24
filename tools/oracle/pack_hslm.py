"""Pack the binary HSLM oracle exports, preserving original Fortran layouts."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exports", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = {"baseline": json.loads(Path("tools/oracle/baseline.json").read_text()), "files": {}}
    for name in [
        "tools/oracle/characterize_hslm.R",
        "tools/oracle/export.R",
        __file__,
        "excavator2/lib/R/LibraryJSLMIn.R",
        "excavator2/lib/F77/FastJointSLMLibraryI.f",
    ]:
        manifest["files"][name] = hashlib.sha256(Path(name).read_bytes()).hexdigest()
    for case in sorted(args.exports.glob("*.RData")):
        arrays = {}
        for folder in case.iterdir():
            metadata = json.loads((folder / "metadata.json").read_text())
            dtype = "<f8" if metadata["type"] == "double" else "<i4"
            values = np.fromfile(folder / "data.bin", dtype=dtype)
            if np.fromfile(folder / "missing.bin", dtype="<i4").any():
                raise ValueError(f"Unexpected missing values: {folder}")
            if metadata["dim"] is not None:
                values = values.reshape(metadata["dim"], order="F")
            arrays[folder.name] = values
        destination = args.output / f"{case.stem}.npz"
        np.savez_compressed(destination, **arrays)
        manifest["files"][str(destination)] = hashlib.sha256(destination.read_bytes()).hexdigest()
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
