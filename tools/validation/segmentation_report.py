"""Compare segmentation on frozen non-normal profiles with fixed FastCall policy."""

import json
from pathlib import Path

import numpy as np

from excavator2.analyze import DEFAULTS, segment_profile, summarize
from excavator2.artifacts import load_manifest
from excavator2.fastcall import assign_labels, fit_fastcall

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests/fixtures/legacy-analysis"


def main():
    target = load_manifest(FIXTURE / "inputs/target", "target")
    prepared = load_manifest(FIXTURE / "inputs/prepared", "prepared")
    report = {
        "comparison": "Frozen legacy HSLM versus current segmentation; same current FastCall",
        "profiles": {},
    }
    for design in ["paired", "pooling"]:
        for sample in ["Test1", "Test2"]:
            old = np.loadtxt(
                FIXTURE / "expected" / design / sample / f"HSLMResults_{sample}.txt",
                dtype=str,
                skiprows=1,
            )
            with np.load(FIXTURE / "inputs/prepared" / prepared["samples"][sample]) as data:
                matrix = data["matrix"]
            rows, _, indices, ids = segment_profile(
                matrix, old[:, 4].astype(float), target, DEFAULTS
            )
            assert np.array_equal(indices, np.arange(len(old)))
            old_values = old[:, 5].astype(float)
            old_starts = np.r_[0, np.flatnonzero(np.diff(old_values) != 0) + 1]
            old_ends = np.r_[old_starts[1:], len(old)]
            starts, ends, values = summarize(rows, ids)
            old_labels = assign_labels(fit_fastcall(old_values[old_starts]).posterior).labels
            labels = assign_labels(fit_fastcall(values).posterior).labels
            old_windows = np.repeat(old_labels, old_ends - old_starts)
            windows = np.repeat(labels, ends - starts)
            changed = np.flatnonzero(windows != old_windows)
            report["profiles"][f"{design}/{sample}"] = {
                "windows": len(old),
                "old_segments": len(old_starts),
                "new_segments": len(starts),
                "old_calls": int(np.count_nonzero(old_labels)),
                "new_calls": int(np.count_nonzero(labels)),
                "changed_window_calls": len(changed),
                "changed_window_indices": changed.tolist(),
                "old_starts": old_starts.tolist(),
                "new_starts": starts.tolist(),
                "old_levels": old_values[old_starts].tolist(),
                "new_levels": values.tolist(),
                "old_labels": old_labels.tolist(),
                "new_labels": labels.tolist(),
            }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
