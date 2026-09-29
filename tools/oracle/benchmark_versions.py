"""Benchmark comparable, plot-free stages after full scientific concordance."""

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
from pathlib import Path

from compare_analysis import compare
from oracle import BASELINE, ROOT, container_command, file_record, save_json
from time_command import measure


def run(folder, repeats):
    folder = folder.resolve()
    evidence = json.loads((folder / "concordance.json").read_text())
    if evidence.get("passed") is not True or evidence["baseline"] != BASELINE:
        raise ValueError("a successful pinned full-pipeline comparison is required")
    if repeats < 3:
        raise ValueError("at least three measured analysis repetitions are required")
    destination = folder / "benchmark"
    destination.mkdir()
    # One scientific analysis worker, with native-library pools also bounded.
    limits = {key: "1" for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]}
    environment = {**os.environ, **limits, "LC_ALL": "C", "TZ": "UTC"}
    samples = []
    for index in range(repeats + 1):  # first pair is an excluded warm-up
        original = destination / f"original-{index}"
        current = destination / f"port-{index}"
        (original / "Results/Test1").mkdir(parents=True)
        old_output = f"/work/benchmark/original-{index}"
        program = "/work/source/excavator2"
        old = container_command(folder) + [
            "env",
            *(f"{key}={value}" for key, value in limits.items()),
            "LC_ALL=C",
            "TZ=UTC",
            "python",
            "/repo/tools/oracle/time_command.py",
            "--report",
            f"/work/benchmark/original-{index}-time.json",
            "--",
            "Rscript",
            f"{program}/lib/R/EXCAVATORInferenceExome.R",
            old_output,
            "/work/legacy-target/synthetic/panel/w_100",
            "/work/paired.yaml",
            "paired",
            program,
            "synthetic",
            "/work/legacy-prepared",
            f"{program}/parameters.yaml",
            "/work/inputs/centromeres.tsv",
        ]
        new = [
            sys.executable,
            "-I",
            "-m",
            "excavator2",
            "analyze",
            "--samples",
            str(folder / "paired.yaml"),
            "--target",
            str(folder / "current-target"),
            "--input",
            str(folder / "current-prepared"),
            "--output",
            str(current),
            "--experiment",
            "paired",
            "--parameters",
            str(folder / "source/excavator2/parameters.yaml"),
        ]
        times = {}
        # Alternate order to avoid always giving the second implementation warmer caches.
        for kind in ["original", "port"] if index % 2 == 0 else ["port", "original"]:
            with (destination / f"{kind}-{index}.log").open("w") as log:
                if kind == "original":
                    subprocess.run(old, stdout=log, stderr=subprocess.STDOUT, check=True)
                    times[kind] = json.loads(
                        (destination / f"original-{index}-time.json").read_text()
                    )["seconds"]
                else:
                    times[kind] = measure(
                        new, cwd=folder, env=environment, stdout=log, stderr=subprocess.STDOUT
                    )
        comparison = compare(original, current)
        if len(comparison["files"]) != 4:
            raise ValueError("benchmark requires all four Test1 scientific outputs")
        if index:
            samples.append({**times, "parity": comparison})
    report = {
        "passed": True,
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "run_url": (
            f"{os.environ.get('GITHUB_SERVER_URL', 'https://github.com')}/"
            f"{os.environ.get('GITHUB_REPOSITORY', 'danilotat/excavator2-py')}/"
            f"actions/runs/{os.environ.get('GITHUB_RUN_ID', '0')}"
        ),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": sys.version,
        "runner": os.environ.get("RUNNER_NAME", "local"),
        "baseline": BASELINE,
        "windows": evidence["prepared"]["Test1"]["windows"],
        "samples": 3,
        "analysis_test_samples": 1,
        "threads": 1,
        "scope": "synthetic target and paired analysis; not full-pipeline or real-WES speedup",
        "exclusions": [
            "plots",
            "preparation timing",
            "Docker startup/image pull",
            "fixture generation",
            "export/conversion/comparison",
        ],
        "timing_policy": (
            "CLI/interpreter startup and artifact I/O included; no cache eviction; "
            "analysis has one excluded warm-up and alternating order"
        ),
        "stages": {
            "target": evidence["timings"]["target"],
            "paired_analysis": {
                "original_seconds": statistics.median(t["original"] for t in samples),
                "port_seconds": statistics.median(t["port"] for t in samples),
                "runs": repeats,
            },
        },
        "analysis_runs": samples,
        "inputs": evidence["inputs"],
        "harness": {
            name: file_record(Path(__file__).with_name(name))
            for name in [
                "benchmark_versions.py",
                "time_command.py",
                "end_to_end.py",
                "run_end_to_end.sh",
            ]
        },
    }
    save_json(destination / "report.json", report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--concordance", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    run(args.concordance, args.repeats)
