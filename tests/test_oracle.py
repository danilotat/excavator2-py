"""Guard the baseline comparator against hiding changed calls or failed runs."""

import importlib.util
import json
from pathlib import Path

import pytest

MODULE = Path(__file__).resolve().parents[1] / "tools/oracle/oracle.py"
spec = importlib.util.spec_from_file_location("oracle", MODULE)
oracle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oracle)


def make_run(root, date="20260924", position=42, status="complete"):
    folder = root / "results/Results/Test1"
    folder.mkdir(parents=True)
    (root / "manifest.json").write_text(
        json.dumps(
            {"status": status, "baseline": {}, "inputs": {}, "source_files": {}, "harness": {}}
        )
    )
    (folder / "calls.vcf").write_text(f"##fileDate={date}\n#CHROM\tPOS\nchr1\t{position}\n")


def test_comparison_ignores_only_date(tmp_path):
    left, right = tmp_path / "left", tmp_path / "right"
    make_run(left)
    make_run(right, date="20260925")
    report = tmp_path / "report.json"
    assert oracle.compare(left, right, report) == 0
    (right / "results/Results/Test1/calls.vcf").write_text("chr1\t43\n")
    assert oracle.compare(left, right, report) == 1


def test_comparison_rejects_failed_run(tmp_path):
    left, right = tmp_path / "left", tmp_path / "right"
    make_run(left)
    make_run(right, status="failed")
    with pytest.raises(ValueError, match="incomplete"):
        oracle.compare(left, right, tmp_path / "report.json")


def test_comparison_rejects_different_inputs(tmp_path):
    left, right = tmp_path / "left", tmp_path / "right"
    make_run(left)
    make_run(right)
    manifest = right / "manifest.json"
    data = json.loads(manifest.read_text())
    data["inputs"] = {"sample.bam": {"sha256": "different"}}
    manifest.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="Different comparison provenance: inputs"):
        oracle.compare(left, right, tmp_path / "report.json")


def test_comparison_detects_missing_artifact(tmp_path):
    left, right = tmp_path / "left", tmp_path / "right"
    make_run(left)
    make_run(right)
    (left / "results/Results/Test1/extra.txt").write_text("chr1\t42\n")
    report = tmp_path / "report.json"
    assert oracle.compare(left, right, report) == 1
    assert json.loads(report.read_text())["differences"] == ["results/Results/Test1/extra.txt"]
