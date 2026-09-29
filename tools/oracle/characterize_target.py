"""M5.1: capture legacy target geometry only; no target-port implementation.

Run with --output pointing to a fresh directory. A pinned source snapshot is run
in Docker; expected successful and failing inputs are recorded, then packed into
small NumPy fixtures. No FASTA/BigWig downloads or feature extraction are involved.
"""

import argparse
import json
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

import numpy as np
from oracle import BASELINE, ROOT, container_command, file_record, save_json

from excavator2.artifacts import read_export


def write_inputs(output):
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    base = [["chr1", 500, 560], ["chr1", 1300, 1400], ["chrY", 800, 850]]
    cases = {
        "standard": (base, True),
        "bare_names": ([[str(c).removeprefix("chr"), s, e] for c, s, e in base], True),
        "unsorted_overlap": (
            [["chr1", 1300, 1400], ["chr1", 500, 560], ["chr1", 550, 700], ["chr1", 550, 650]],
            True,
        ),
        "gap_endpoints": (
            [
                ["chr1", 500, 560],
                ["chr1", 1200, 1700],
                ["chr1", 1400, 1500],
                ["chr1", 1500, 1550],
                ["chr1", 1600, 1650],
                ["chr1", 2000, 2100],
            ],
            True,
        ),
        "coordinate_names_ignored": (base, True),
        "no_eligible_gap": ([["chr1", 10, 2990]], False),
        "no_alternate_gap": (base, False),
        "missing_chromosome_gap": (base, False),
    }
    for name, (bed, _) in cases.items():
        folder = output / "cases" / name
        folder.mkdir(parents=True)
        coords = [
            [c.removeprefix("chr") if name == "bare_names" else c, 0, 3000] for c in chromosomes
        ]
        if name == "coordinate_names_ignored":
            for row, replacement in zip(coords, reversed(chromosomes), strict=True):
                row[0] = replacement
        gaps = [[i, c, 2400, 2450] for i, c in enumerate(chromosomes)]
        if name == "gap_endpoints":
            gaps[0][2:] = [1500, 1600]
        if name == "missing_chromosome_gap":
            gaps = [row for row in gaps if row[1] != "chr2"]
        if name != "no_alternate_gap":
            gaps.append([99, "chr1_alt", 10, 20])
        for filename, rows in [
            ("target.bed", bed),
            ("chromosomes.tsv", coords),
            ("gaps.tsv", gaps),
        ]:
            header = "id\tchrom\tstart\tend\n" if filename == "gaps.tsv" else ""
            (folder / filename).write_text(
                header + "".join("\t".join(map(str, row)) + "\n" for row in rows)
            )
    save_json(
        output / "cases.json", {name: {"expected_success": ok} for name, (_, ok) in cases.items()}
    )
    return cases


def run(output, compare=False):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    cases = write_inputs(output)
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
    with (output / "container.log").open("w") as log:
        subprocess.run(
            container_command(output)
            + ["Rscript", "/repo/tools/oracle/characterize_target.R", "/work"],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
        subprocess.run(
            container_command(output)
            + ["Rscript", "/repo/tools/oracle/export.R", "/work/cases", "/work/exports"],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    results = json.loads((output / "results.json").read_text())
    fixtures = output / "fixtures"
    fixtures.mkdir()
    manifest = {
        "baseline": BASELINE,
        "scope": "FilterTarget.R geometry only; window=100, flank=200; feature extraction deferred",
        "source": file_record(source / "excavator2/lib/R/FilterTarget.R"),
        "harness": {
            name: file_record(Path(__file__).with_name(name))
            for name in ["characterize_target.py", "characterize_target.R", "export.R"]
        },
        "cases": {},
    }
    for name, (_, expected_success) in cases.items():
        folder = output / "cases" / name
        status = results[name]["status"]
        if (status == 0) != expected_success:
            raise ValueError(
                f"unexpected legacy result for {name}: status {status}; see {folder / 'stderr.log'}"
            )
        destination = fixtures / name
        destination.mkdir()
        for filename in ["target.bed", "chromosomes.tsv", "gaps.tsv", "stderr.log"]:
            shutil.copyfile(folder / filename, destination / filename)
        record = {"status": status, "expected_success": expected_success, "files": {}}
        if expected_success:
            target = read_export(output / "exports" / name / "output/panel.RData/MyTarget")
            text_target = np.loadtxt(
                folder / "output/Filtered.txt", dtype=str, delimiter="\t", ndmin=2
            )
            np.testing.assert_array_equal(target, text_target)
            names = (folder / "output/panel_chromosome.txt").read_text().split()
            np.savez_compressed(
                destination / "expected.npz", target=target, chromosomes=np.array(names)
            )
            record["rows"] = len(target)
        for path in destination.iterdir():
            record["files"][path.name] = file_record(path)
        manifest["cases"][name] = record
    save_json(fixtures / "manifest.json", manifest)
    if compare:
        compare_current(fixtures)
    print(
        json.dumps(
            {
                name: {k: v for k, v in record.items() if k != "files"}
                for name, record in manifest["cases"].items()
            },
            indent=2,
        )
    )


def compare_current(fixtures):
    from excavator2.target import target_geometry

    manifest = json.loads((fixtures / "manifest.json").read_text())
    for name, record in manifest["cases"].items():
        folder = fixtures / name
        try:
            actual = target_geometry(
                folder / "target.bed", folder / "chromosomes.tsv", folder / "gaps.tsv", 100
            )
        except ValueError:
            if record["expected_success"]:
                raise
            continue
        if not record["expected_success"]:
            raise AssertionError(f"port accepted legacy geometry failure {name}")
        with np.load(folder / "expected.npz", allow_pickle=False) as expected:
            np.testing.assert_array_equal(actual, expected["target"], err_msg=name)
            np.testing.assert_array_equal(
                list(dict.fromkeys(actual[:, 0])), expected["chromosomes"]
            )
    save_json(
        fixtures.parent / "concordance.json",
        {
            "passed": True,
            "cases": len(manifest["cases"]),
            "comparison": "exact geometry and success/failure outcomes",
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compare-current", action="store_true")
    args = parser.parse_args()
    run(args.output, args.compare_current)


if __name__ == "__main__":
    main()
