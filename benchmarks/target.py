"""Fresh-process target timings and exact artifact comparison against a baseline.

Run with each installed Python under comparison. Input caches are not evicted;
these are repeated/warm-cache measurements, not cold-storage benchmarks.
"""

import argparse
import hashlib
import json
import platform
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from excavator2 import artifacts, target, target_features, target_pipeline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-target", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    if args.worker:
        wall, cpu = time.perf_counter(), time.process_time()
        target_pipeline.run_target(args.config, args.output / "target")
        elapsed, cpu_elapsed = time.perf_counter() - wall, time.process_time() - cpu
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # macOS reports bytes; Linux reports KiB. Capture before comparator work.
        peak_bytes = rss if sys.platform == "darwin" else rss * 1024
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools/oracle"))
        from compare_target import compare

        comparison = compare(args.expected_target, args.output / "target")
        manifest = artifacts.load_manifest(args.output / "target", "target")
        with np.load(args.output / "target" / manifest["preparation"]) as data:
            windows = len(data["target"])
        result = {
            "wall_seconds": elapsed,
            "cpu_seconds": cpu_elapsed,
            "peak_rss_bytes": peak_bytes,
            "windows": windows,
            "windows_per_second": windows / elapsed,
            "exact_target_parity": comparison["passed"],
            "input_files": manifest["inputs"],
        }
    else:
        runs = []
        for index in range(args.repeats):
            folder = args.output / f"run-{index}"
            subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--worker",
                    "--config",
                    str(args.config.resolve()),
                    "--output",
                    str(folder.resolve()),
                    "--expected-target",
                    str(args.expected_target.resolve()),
                ],
                check=True,
            )
            runs.append(json.loads((folder / "report.json").read_text()))
        result = {
            "platform": platform.platform(),
            "python": sys.version,
            "numpy": np.__version__,
            "threads": 1,
            "cache_policy": "No cache eviction; repeated local-file runs, not cold I/O",
            "timing_scope": "run_target including reference reads, hashing and artifact writes",
            "source_sha256": {
                Path(m.__file__).name: artifacts.digest(m.__file__)
                for m in [artifacts, target, target_features, target_pipeline]
            },
            "benchmark_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "runs": runs,
            "median_wall_seconds": statistics.median(r["wall_seconds"] for r in runs),
            "median_cpu_seconds": statistics.median(r["cpu_seconds"] for r in runs),
            "median_peak_rss_bytes": statistics.median(r["peak_rss_bytes"] for r in runs),
        }
    (args.output / "report.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
