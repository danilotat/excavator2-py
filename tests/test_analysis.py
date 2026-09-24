"""Exercise the installed CLI against non-normal original-R paired/pooling outputs."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from excavator2.analyze import experimental_design, run_analysis, summarize
from excavator2.artifacts import load_manifest, read_yaml

FIXTURE = Path(__file__).parent / "fixtures" / "legacy-analysis"


def lines(path):
    return [line for line in path.read_text().splitlines() if not line.startswith("##fileDate=")]


@pytest.mark.parametrize("experiment", ["paired", "pooling"])
def test_cli_matches_original_non_normal_outputs(tmp_path, experiment):
    output = tmp_path / "results"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "excavator2",
            "analyze",
            "--samples",
            str(FIXTURE / "samples.yaml"),
            "--input",
            str(FIXTURE / "inputs/prepared"),
            "--target",
            str(FIXTURE / "inputs/target"),
            "--output",
            str(output),
            "--experiment",
            experiment,
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    for sample in ["Test1", "Test2"]:
        for expected in (FIXTURE / "expected" / experiment / sample).iterdir():
            assert lines(output / "Results" / sample / expected.name) == lines(expected)
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["samples"]["Test1"]["calls"] == 4


def test_failed_analysis_does_not_publish_partial_results(tmp_path):
    parameters = tmp_path / "parameters.yaml"
    parameters.write_text(
        "HSLM: {Omega: 0, Theta: 0.1, D_norm: 100}\n"
        "FastCall: {Cellularity: 1, d: 0.5, u: 0.35, minExons: 4}\n"
    )
    output = tmp_path / "results"
    with pytest.raises(ValueError, match="omega"):
        run_analysis(
            FIXTURE / "samples.yaml",
            FIXTURE / "inputs/prepared",
            FIXTURE / "inputs/target",
            output,
            "paired",
            parameters,
        )
    assert not output.exists()
    assert not list(tmp_path.glob(".excavator2-*"))


def test_manifest_checksums_and_target_identity(tmp_path):
    prepared = tmp_path / "prepared"
    shutil.copytree(FIXTURE / "inputs/prepared", prepared)
    manifest = json.loads((prepared / "manifest.json").read_text())
    manifest["target_id"] = "wrong-target"
    (prepared / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="identities"):
        run_analysis(
            FIXTURE / "samples.yaml",
            prepared,
            FIXTURE / "inputs/target",
            tmp_path / "output",
            "paired",
        )
    (prepared / next(iter(manifest["files"]))).write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="checksum"):
        load_manifest(prepared, "prepared")


def test_design_order_and_duplicate_yaml(tmp_path):
    samples = {"T2": "test2", "C2": "control2", "C1": "control1", "T1": "test1"}
    assert experimental_design(samples, "paired") == [
        ("test1", ["control1"]),
        ("test2", ["control2"]),
    ]
    assert experimental_design(samples, "pooled") == [
        ("test2", ["control2", "control1"]),
        ("test1", ["control2", "control1"]),
    ]
    with pytest.raises(ValueError, match="matching"):
        experimental_design({"C2": "control", "T1": "test"}, "paired")
    path = tmp_path / "samples.yaml"
    path.write_text("T1: one\nT1: two\n")
    with pytest.raises(ValueError, match="duplicate"):
        read_yaml(path)


def test_makedata_preserves_cross_chromosome_grouping():
    rows = np.array(
        [["chr1", "1", "1", "2", "0", "0", "IN"], ["chr2", "3", "3", "4", "0", "0", "OUT"]]
    )
    starts, ends, values = summarize(rows)
    assert_array_equal(starts, [0])
    assert_array_equal(ends, [2])
    assert_array_equal(values, [0])


def test_legacy_converter_preserves_character_values(tmp_path):
    from excavator2.artifacts import convert_legacy

    def export(matrix, folder):
        folder.mkdir(parents=True)
        (folder / "metadata.json").write_text(
            json.dumps(
                {
                    "type": "character",
                    "length": matrix.size,
                    "dim": list(matrix.shape),
                    "order": "F",
                    "endian": "little",
                }
            )
        )
        (folder / "values.json").write_text(json.dumps(matrix.ravel(order="F").tolist()))
        np.zeros(matrix.size, dtype="<i4").tofile(folder / "missing.bin")

    with np.load(FIXTURE / "inputs/prepared/sample-0.npz") as f:
        matrix = f["matrix"]
    export(matrix, tmp_path / "exports/prepared/Control/RCNorm/Control.NRC.RData/MatrixNorm")
    for index, chromosome in enumerate(["chr1", "chr2"]):
        with np.load(FIXTURE / f"inputs/target/reference-{index}.npz") as f:
            export(f["matrix"], tmp_path / f"exports/target/FRB/FRB.{chromosome}.RData/FRBData")
    (tmp_path / "chromosomes.txt").write_text("chr1 chr2")
    (tmp_path / "centromeres.txt").write_text(
        "chrom\tstart\tend\nchr1\t7000\t11000\nchr2\t7000\t11000\n"
    )
    output = tmp_path / "converted"
    convert_legacy(
        tmp_path / "exports/prepared",
        tmp_path / "exports/target",
        tmp_path / "chromosomes.txt",
        tmp_path / "centromeres.txt",
        "synthetic",
        output,
    )
    prepared = load_manifest(output / "prepared", "prepared")
    target = load_manifest(output / "target", "target")
    assert prepared["target_id"] == target["target_id"]
    with np.load(output / "prepared/sample-0.npz") as f:
        assert_array_equal(f["matrix"], matrix)
