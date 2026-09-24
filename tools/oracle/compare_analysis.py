"""Compare analysis outputs with exact decisions and bounded continuous drift."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def compare(original, rewritten):
    reports = {}
    for expected in sorted(original.glob("Results/*/*")):
        if expected.suffix not in (".txt", ".vcf"):
            continue
        relative = expected.relative_to(original)
        actual = rewritten / relative
        old = expected.read_text().splitlines()
        new = actual.read_text().splitlines()
        if expected.suffix == ".vcf":
            old = [s for s in old if not s.startswith("##fileDate=")]
            new = [s for s in new if not s.startswith("##fileDate=")]
        report = {
            "original_sha256": hashlib.sha256(expected.read_bytes()).hexdigest(),
            "rewritten_sha256": hashlib.sha256(actual.read_bytes()).hexdigest(),
            "rows": len(old),
            "text_equal": old == new,
        }
        if expected.name.startswith("HSLMResults"):
            if old[0] != new[0] or len(old) != len(new):
                raise ValueError(f"HSLM header/row count differs: {relative}")
            a = np.array([s.split("\t") for s in old[1:]])
            b = np.array([s.split("\t") for s in new[1:]])
            if not np.array_equal(a[:, [0, 6]], b[:, [0, 6]]) or not np.array_equal(
                a[:, 1:4].astype(float), b[:, 1:4].astype(float)
            ):
                raise ValueError(f"HSLM metadata differs: {relative}")
            x, y = a[:, 4:6].astype(float), b[:, 4:6].astype(float)
            if not np.allclose(x, y, rtol=1e-13, atol=2e-15):
                raise ValueError(f"HSLM numeric drift: {relative}")
            if not np.array_equal(np.diff(x[:, 1]) != 0, np.diff(y[:, 1]) != 0):
                raise ValueError(f"segment boundaries differ: {relative}")
            report["max_absolute_error"] = float(np.max(np.abs(x - y)))
            report["decision_mismatches"] = 0
        elif old != new:
            raise ValueError(f"call table/VCF differs: {relative}")
        reports[str(relative)] = report
    if not reports:
        raise ValueError("no expected outputs found")
    return {"files": reports, "excluded_fields": ["VCF fileDate"], "passed": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path)
    parser.add_argument("rewritten", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    args.report.write_text(json.dumps(compare(args.original, args.rewritten), indent=2) + "\n")


if __name__ == "__main__":
    main()
