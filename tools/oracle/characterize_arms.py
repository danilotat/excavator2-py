"""Capture original full-inference outcomes at chromosome-arm boundaries."""

import argparse
import subprocess
import tarfile
import tempfile
from pathlib import Path

from oracle import BASELINE, ROOT, container_command, file_record, save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compare-fixtures", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
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
    with (output / "legacy.log").open("w") as log:
        subprocess.run(
            container_command(output) + ["bash", "/repo/tools/oracle/run_concordance.sh"],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
        subprocess.run(
            container_command(output)
            + [
                "Rscript",
                "/repo/tools/oracle/characterize_arms.R",
                "/work/legacy",
                "/work/source/excavator2",
                "/work/arms",
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    # Keep only boundary metadata, original diagnostics and HSLM output evidence.
    evidence = output / "evidence"
    evidence.mkdir()
    for path in sorted((output / "arms").rglob("*")):
        if not path.is_file() or not (
            path.name in ("cases.tsv", "legacy.log") or path.name.startswith("HSLMResults_")
        ):
            continue
        destination = evidence / path.relative_to(output / "arms")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(path.read_bytes())
    files = {
        str(p.relative_to(evidence)): file_record(p) for p in evidence.rglob("*") if p.is_file()
    }
    if args.compare_fixtures:
        fixtures = ROOT / "tests/fixtures/legacy-analysis/arms"
        expected = {
            str(p.relative_to(fixtures))
            for p in fixtures.rglob("*")
            if p.is_file() and p.name != "manifest.json"
        }
        if set(files) != expected:
            raise ValueError("arm evidence file inventory differs")
        for name in files:
            if (evidence / name).read_bytes() != (fixtures / name).read_bytes():
                raise ValueError(f"arm evidence differs: {name}")
    save_json(
        evidence / "manifest.json",
        {
            "baseline": BASELINE,
            "source": {
                name: file_record(source / "excavator2" / name)
                for name in ["lib/R/EXCAVATORInferenceExome.R", "lib/R/LibraryJSLMIn.R"]
            },
            "harness": {
                name: file_record(ROOT / "tools/oracle" / name)
                for name in ["characterize_arms.R", "characterize_analysis.R", "run_concordance.sh"]
            },
            "files": files,
        },
    )


if __name__ == "__main__":
    main()
