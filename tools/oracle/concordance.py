"""Run old and current analysis on fresh shared inputs, then compare both designs.

Use an installed-wheel Python to qualify that wheel. Docker is required, but no
BAM/FASTA downloads or existing .oracle state are used. Any failure exits nonzero
and leaves logs, inputs, outputs and a machine-readable report in --output.
"""

import argparse
import json
import platform
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

from compare_analysis import compare
from compare_preparation import compare_preparation
from make_preparation_bams import make_bams
from oracle import BASELINE, ROOT, container_command, file_record, save_json

from excavator2 import __version__, _core
from excavator2.artifacts import convert_legacy


def run(output):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = {
        "passed": False,
        "baseline": BASELINE,
        "current_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "python": sys.version,
        "platform": platform.platform(),
        "package_version": __version__,
        "extension": file_record(Path(_core.__file__)),
        "harness": {
            name: file_record(Path(__file__).with_name(name))
            for name in [
                "concordance.py",
                "run_concordance.sh",
                "characterize_analysis.R",
                "export.R",
                "compare_analysis.py",
                "compare_preparation.py",
                "make_preparation_bams.py",
                "characterize_preparation_pipeline.R",
                "run_preparation_concordance.sh",
            ]
        },
        "designs": {},
    }
    try:
        source = output / "source"
        source.mkdir()
        with tempfile.TemporaryFile() as archive:
            subprocess.run(
                ["git", "archive", BASELINE["source_commit"], "excavator2"],
                cwd=ROOT,
                stdout=archive,
                check=True,
            )
            archive.seek(0)
            with tarfile.open(fileobj=archive) as tar:
                tar.extractall(source, filter="data")
        report["source_files"] = {
            str(path.relative_to(source)): file_record(path)
            for path in sorted(source.rglob("*"))
            if path.is_file()
        }
        with (output / "legacy.log").open("w") as log:
            subprocess.run(
                container_command(output) + ["bash", "/repo/tools/oracle/run_concordance.sh"],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        convert_legacy(
            output / "exports/prepared",
            output / "exports/target",
            output / "legacy/target/synthetic_chromosome.txt",
            output / "legacy/centromeres.txt",
            "synthetic",
            output / "converted",
        )
        for mode in ["paired", "pooling"]:
            try:
                with (output / f"current-{mode}.log").open("w") as log:
                    subprocess.run(
                        [
                            sys.executable,
                            "-m",
                            "excavator2",
                            "analyze",
                            "--samples",
                            str(output / "legacy/paired.yaml"),
                            "--input",
                            str(output / "converted/prepared"),
                            "--target",
                            str(output / "converted/target"),
                            "--output",
                            str(output / f"current-{mode}"),
                            "--experiment",
                            mode,
                        ],
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        check=True,
                    )
                result = compare(output / "legacy" / mode, output / f"current-{mode}")
                if len(result["files"]) != 8:
                    raise ValueError("expected four output files for each of two samples")
                report["designs"][mode] = result
            except (ValueError, OSError, subprocess.CalledProcessError) as error:
                report["designs"][mode] = {"passed": False, "error": str(error)}
        make_bams(output)
        with (output / "legacy-preparation.log").open("w") as log:
            subprocess.run(
                container_command(output)
                + ["bash", "/repo/tools/oracle/run_preparation_concordance.sh"],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        convert_legacy(
            None,
            output / "exports/preparation-target",
            output / "legacy/target/synthetic_chromosome.txt",
            output / "legacy/centromeres.txt",
            "synthetic",
            output / "preparation-input",
        )
        with (output / "current-preparation.log").open("w") as log:
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "excavator2",
                    "prepare",
                    "--samples",
                    str(output / "bam-samples.yaml"),
                    "--target",
                    str(output / "preparation-input/target"),
                    "--output",
                    str(output / "current-prepared"),
                    "--threads",
                    "2",
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        report["preparation"] = compare_preparation(
            output / "exports/legacy-prepared", output / "current-prepared"
        )
        for mode in ["paired", "pooling"]:
            with (output / f"current-prepared-{mode}.log").open("w") as log:
                subprocess.run(
                    [
                        sys.executable,
                        "-m",
                        "excavator2",
                        "analyze",
                        "--samples",
                        str(output / "legacy/paired.yaml"),
                        "--input",
                        str(output / "current-prepared"),
                        "--target",
                        str(output / "preparation-input/target"),
                        "--output",
                        str(output / f"current-prepared-{mode}"),
                        "--experiment",
                        mode,
                    ],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            report["designs"][f"prepared-{mode}"] = compare(
                output / f"legacy-prepared-{mode}", output / f"current-prepared-{mode}"
            )
        report["passed"] = all(result["passed"] for result in report["designs"].values())
    except (ValueError, OSError, subprocess.CalledProcessError, tarfile.TarError) as error:
        report["error"] = str(error)
    finally:
        save_json(output / "concordance.json", report)
    print(json.dumps({"passed": report["passed"], "report": str(output / "concordance.json")}))
    return 0 if report["passed"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.output))


if __name__ == "__main__":
    main()
