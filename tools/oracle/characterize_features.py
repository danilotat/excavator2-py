"""Capture original TargetCreate.sh features on tiny FASTA/BigWig references."""

import argparse
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pyBigWig
import pysam
from oracle import BASELINE, ROOT, container_command, file_record, save_json

from excavator2.artifacts import read_export


def run(output, compare=False):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    reference = output / "reference.fa"
    reference.write_text(">chr1\n" + "ACGTacgtNNRYCGATNnGC" * 5 + "\n>chr2\n" + "N" * 100 + "\n")
    pysam.faidx(str(reference))
    with pyBigWig.open(str(output / "reference.bw"), "w") as bw:
        bw.addHeader([("chr1", 100), ("chr2", 100)])
        bw.addEntries(
            ["chr1"] * 4, [0, 5, 20, 70], ends=[3, 10, 40, 100], values=[0.1, 0.7, 0.0, 0.123456789]
        )
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
    intervals = [
        ("chr1", 1, 3),
        ("chr1", 2, 9),
        ("chr1", 5, 7),
        ("chr1", 9, 12),
        ("chr1", 9, 19),
        ("chr1", 21, 30),
        ("chr1", 41, 50),
        ("chr1", 71, 99),
        ("chr2", 1, 99),
    ]
    cases = {
        "standard": intervals,
        "bare_names": [(c.removeprefix("chr"), s, e) for c, s, e in intervals],
        "past_end": [("chr1", 95, 105)],
        "last_base": [("chr1", 100, 100)],
        "partial_skips": [("chr1", 1, 3), ("chr1", 100, 100), ("chr3", 1, 9)],
        "zero_start": [("chr1", 0, 9)],
        "missing_contig": [("chr3", 1, 9)],
    }
    records = {}
    for name, rows in cases.items():
        folder = output / name
        folder.mkdir()
        for sub in ["GCC", "MAP", "FRB"]:
            (folder / sub).mkdir()
        target = np.array(
            [(c, str(s), str(e), f"a{i + 1}", "IN") for i, (c, s, e) in enumerate(rows)]
        )
        np.savetxt(folder / "target.tsv", target, fmt="%s", delimiter="\t")
        (folder / "panel_chromosome.txt").write_text(" ".join(dict.fromkeys(target[:, 0])) + "\n")
        with (folder / "legacy.log").open("w") as log:
            result = subprocess.run(
                container_command(output)
                + [
                    "bash",
                    "/work/source/excavator2/lib/bash/TargetCreate.sh",
                    "/work/reference.bw",
                    f"/work/{name}/target.tsv",
                    "/work/source/excavator2",
                    "unused",
                    "unused",
                    "/work/reference.fa",
                    f"/work/{name}",
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        expected_success = name in {"standard", "bare_names", "past_end"}
        if (result.returncode == 0) != expected_success:
            raise ValueError(f"unexpected legacy status for {name}: {result.returncode}")
        records[name] = {"status": result.returncode}
    subprocess.run(
        container_command(output)
        + ["Rscript", "/repo/tools/oracle/export.R", "/work", "/work/exports"],
        check=True,
    )
    fixtures = output / "fixtures"
    fixtures.mkdir()
    for name in ["reference.fa", "reference.fa.fai", "reference.bw"]:
        shutil.copyfile(output / name, fixtures / name)
    for name in cases:
        folder = fixtures / name
        folder.mkdir()
        for filename in ["target.tsv", "legacy.log"]:
            shutil.copyfile(output / name / filename, folder / filename)
        if records[name]["status"] == 0:
            arrays = {}
            target = np.loadtxt(folder / "target.tsv", dtype=str, ndmin=2)
            for chrom in dict.fromkeys(target[:, 0]):
                for kind, obj in [("GCC", "GCContent"), ("MAP", "MapMed"), ("FRB", "FRBData")]:
                    filename = f"{'Map' if kind == 'MAP' else kind}.{chrom}.RData"
                    arrays[f"{kind}_{chrom}"] = read_export(
                        output / "exports" / name / kind / filename / obj
                    )
            np.savez_compressed(folder / "expected.npz", **arrays)
    save_json(
        fixtures / "manifest.json",
        {
            "baseline": BASELINE,
            "cases": records,
            "harness": file_record(Path(__file__)),
            "source": {
                str(p.relative_to(source)): file_record(p)
                for p in [
                    source / "excavator2/lib/bash/TargetCreate.sh",
                    *[source / f"excavator2/lib/R/Save{k}.R" for k in ["Map", "GCC", "FRB"]],
                ]
            },
            "files": {
                str(p.relative_to(fixtures)): file_record(p)
                for p in fixtures.rglob("*")
                if p.is_file()
            },
        },
    )
    if compare:
        compare_current(fixtures)
    print(records)


def compare_current(fixtures):
    """Fail on any feature difference or changed success/failure outcome."""
    import json

    from excavator2.target_features import target_features

    manifest = json.loads((fixtures / "manifest.json").read_text())
    for name, record in manifest["cases"].items():
        target = np.loadtxt(fixtures / name / "target.tsv", dtype=str, ndmin=2)
        try:
            actual = target_features(target, fixtures / "reference.fa", fixtures / "reference.bw")
        except ValueError:
            if record["status"] == 0:
                raise
            continue
        if record["status"] != 0:
            raise AssertionError(f"port accepted legacy failure {name}")
        with np.load(fixtures / name / "expected.npz", allow_pickle=False) as expected:
            keys = {f"{kind}_{c}" for c in actual for kind in ["GCC", "MAP", "FRB"]}
            if keys != set(expected.files):
                raise AssertionError(f"different feature objects for {name}")
            for chromosome, values in actual.items():
                for kind, key in [("GCC", "gc"), ("MAP", "mappability"), ("FRB", "first_base")]:
                    np.testing.assert_array_equal(
                        values[key],
                        expected[f"{kind}_{chromosome}"],
                        err_msg=f"{name}/{chromosome}/{kind}",
                    )
    save_json(
        fixtures.parent / "concordance.json",
        {
            "status": "passed",
            "cases": len(manifest["cases"]),
            "comparison": "exact feature arrays and success/failure outcomes",
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--compare-current", action="store_true")
    args = parser.parse_args()
    run(args.output, args.compare_current)
