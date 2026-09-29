"""Compare a new target artifact with a converted original target, exactly."""

import argparse
import json
from pathlib import Path

import numpy as np

from excavator2.artifacts import load_manifest


def compare(original, rewritten):
    original, rewritten = Path(original), Path(rewritten)
    old, new = load_manifest(original, "target"), load_manifest(rewritten, "target")
    for key in ["target_id", "assembly", "chromosomes", "centromeres"]:
        if old[key] != new[key]:
            raise ValueError(f"target manifest differs: {key}")
    reports = {}
    pairs = [("preparation", old["preparation"], new["preparation"])]
    pairs.extend(
        (chrom, old["references"][chrom], new["references"][chrom]) for chrom in old["chromosomes"]
    )
    for label, old_file, new_file in pairs:
        with (
            np.load(original / old_file, allow_pickle=False) as x,
            np.load(rewritten / new_file, allow_pickle=False) as y,
        ):
            if set(x.files) != set(y.files):
                raise ValueError(f"different array keys: {label}")
            for key in x.files:
                np.testing.assert_array_equal(x[key], y[key], err_msg=f"{label}/{key}")
                reports[f"{label}/{key}"] = {"shape": list(x[key].shape), "exact": True}
    return {"passed": True, "target_id": new["target_id"], "arrays": reports}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path)
    parser.add_argument("rewritten", type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    args.report.write_text(json.dumps(compare(args.original, args.rewritten), indent=2) + "\n")
