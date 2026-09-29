"""Compare M4 counts and normalized matrices against the saved supplied-data oracle."""

import argparse
import json
from pathlib import Path

import numpy as np

from excavator2.artifacts import digest, load_manifest, read_export


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("legacy_exports", type=Path)
    parser.add_argument("prepared", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    prepared = load_manifest(args.prepared, "prepared")
    target = load_manifest(args.target, "target")
    report = {
        "baseline": json.loads(Path(__file__).with_name("baseline.json").read_text()),
        "prepared_manifest_sha256": digest(args.prepared / "manifest.json"),
        "samples": {},
    }
    for sample, filename in prepared["samples"].items():
        root = args.legacy_exports / sample
        counts = np.concatenate(
            [
                read_export(root / "RC" / f"{sample}.RC.{c}.RData" / "RC")
                for c in target["chromosomes"]
            ]
        )
        old = read_export(root / "RCNorm" / f"{sample}.NRC.RData" / "MatrixNorm")
        with np.load(args.prepared / filename) as current:
            np.testing.assert_array_equal(counts, current["counts"])
            np.testing.assert_array_equal(
                old[:, [0, 1, 2, 3, 4, 6]], current["matrix"][:, [0, 1, 2, 3, 4, 6]]
            )
            expected, actual = old[:, 5].astype(float), current["matrix"][:, 5].astype(float)
            np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=2e-15)
            report["samples"][sample] = {
                "windows": len(counts),
                "count_mismatches": 0,
                "normalized_max_absolute_error": float(max(abs(actual - expected))),
                "normalized_text_equal": bool(np.array_equal(old[:, 5], current["matrix"][:, 5])),
            }
    report["passed"] = True
    args.report.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
