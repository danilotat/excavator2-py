"""Prove the live-CI comparator rejects changed calls and hidden output loss."""

import importlib.util
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "compare_analysis", ROOT / "tools/oracle/compare_analysis.py"
)
comparator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparator)


@pytest.fixture
def runs(tmp_path):
    expected = ROOT / "tests/fixtures/legacy-analysis/expected/paired"
    original, current = tmp_path / "old", tmp_path / "current"
    for folder in [original, current]:
        shutil.copytree(expected, folder / "Results")
    return original, current


def test_comparison_accepts_only_date_and_bounded_continuous_drift(runs):
    original, current = runs
    path = current / "Results/Test1/EXCAVATORRegionCall_Test1.vcf"
    lines = path.read_text().splitlines()
    lines[1] = "##fileDate=20000101"
    path.write_text("\n".join(lines) + "\n")
    path = current / "Results/Test1/HSLMResults_Test1.txt"
    lines = path.read_text().splitlines()
    fields = lines[1].split("\t")
    fields[4] = str(float(fields[4]) + 1e-16)
    lines[1] = "\t".join(fields)
    path.write_text("\n".join(lines) + "\n")
    assert comparator.compare(original, current)["passed"]


@pytest.mark.parametrize(
    "mutation", ["missing", "extra", "call", "vcf", "coordinate", "boundary", "drift"]
)
def test_comparison_rejects_changed_results(runs, mutation):
    original, current = runs
    folder = current / "Results/Test1"
    if mutation == "missing":
        (folder / "FastCallResults_Test1.txt").unlink()
    elif mutation == "extra":
        (folder / "unexpected.txt").write_text("unexpected")
    elif mutation == "vcf":
        path = folder / "EXCAVATORRegionCall_Test1.vcf"
        path.write_text(path.read_text().replace("1/1:", "0/1:", 1))
    else:
        path = folder / (
            "FastCallResults_Test1.txt" if mutation == "call" else "HSLMResults_Test1.txt"
        )
        lines = path.read_text().splitlines()
        fields = lines[1].split("\t")
        if mutation == "call":
            fields[6] = "0"
        elif mutation == "coordinate":
            fields[2] = str(float(fields[2]) + 1)
        elif mutation == "boundary":
            # Smaller than numerical tolerance, but creates a new segment boundary.
            fields[5] = str(float(fields[5]) + 1e-15)
        else:
            fields[4] = str(float(fields[4]) + 0.1)
        lines[1] = "\t".join(fields)
        path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError):
        comparator.compare(original, current)
