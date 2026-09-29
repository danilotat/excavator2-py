"""Compare exact counts and all normalization checkpoints against fresh R exports."""

import numpy as np

from excavator2.artifacts import load_manifest, read_export


def compare_preparation(original, prepared):
    manifest = load_manifest(prepared, "prepared")
    reports = {}
    for sample, filename in manifest["samples"].items():
        with np.load(prepared / filename) as current:
            trace = original / sample / "trace.RData"
            counts = read_export(trace / "counts")
            if not np.array_equal(counts, current["counts"]):
                raise ValueError(f"{sample}: changed integer read counts")
            report = {"count_mismatches": 0, "windows": len(counts), "normalization": {}}
            for field, legacy in [
                ("weighted", "weighted"),
                ("size", "size"),
                ("map", "mapped"),
                ("gc", "corrected"),
                ("normalized", "normalized"),
            ]:
                old, new = read_export(trace / legacy), current[field]
                if not np.allclose(old, new, rtol=1e-13, atol=2e-15):
                    raise ValueError(f"{sample}: changed {field} correction")
                report["normalization"][field] = float(np.max(np.abs(old - new)))
            old_matrix = read_export(
                original / sample / "RCNorm" / f"{sample}.NRC.RData" / "MatrixNorm"
            )
            if not np.array_equal(
                old_matrix[:, [0, 1, 2, 3, 4, 6]], current["matrix"][:, [0, 1, 2, 3, 4, 6]]
            ):
                raise ValueError(f"{sample}: changed prepared window metadata")
            if not np.allclose(
                old_matrix[:, 5].astype(float),
                current["matrix"][:, 5].astype(float),
                rtol=1e-13,
                atol=2e-15,
            ):
                raise ValueError(f"{sample}: changed prepared normalized counts")
            reports[sample] = report
    if set(reports) != {path.name for path in original.iterdir() if path.is_dir()}:
        raise ValueError("prepared sample inventory differs")
    return {"passed": True, "samples": reports}
