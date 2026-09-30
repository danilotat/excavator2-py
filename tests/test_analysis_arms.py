"""Full legacy inference evidence for successful and unsupported arm geometry."""

import copy
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2.analyze import DEFAULTS, segment_profile
from excavator2.artifacts import load_manifest

FIXTURE = Path(__file__).parent / "fixtures/legacy-analysis"
EDGES = FIXTURE / "arms"
with (EDGES / "cases.tsv").open() as handle:
    CASES = list(csv.DictReader(handle, delimiter="\t"))

ERRORS = {
    "no_long": ("no long arm", "NA/NaN argument"),
    "endpoints": ("excludes windows", "number of items to replace"),
    "inside": ("excludes windows", "number of items to replace"),
}


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["case"])
@pytest.mark.parametrize("sample", ["Test1", "Test2"])
def test_arm_geometry_against_original_inference(case, sample):
    target = copy.deepcopy(load_manifest(FIXTURE / "inputs/target", "target"))
    target["centromeres"]["chr1"] = [int(case["start"]), int(case["end"])]
    prepared = load_manifest(FIXTURE / "inputs/prepared", "prepared")

    def matrix(name):
        with np.load(FIXTURE / "inputs/prepared" / prepared["samples"][name]) as data:
            return data["matrix"]

    test = matrix(sample)
    values = np.loadtxt(
        FIXTURE / "expected/paired" / sample / f"HSLMResults_{sample}.txt",
        dtype=str,
        skiprows=1,
    )[:, 4].astype(float)
    if case["case"] in ERRORS:
        port_error, legacy_error = ERRORS[case["case"]]
        assert int(case["status"]) != 0
        assert legacy_error in (EDGES / case["case"] / "legacy.log").read_text()
        with pytest.raises(ValueError, match=port_error):
            segment_profile(test, values, target, DEFAULTS)
        return
    actual, _, indices, ids = segment_profile(test, values, target, DEFAULTS)
    assert_array_equal(indices, np.arange(len(test)))
    assert_array_equal(actual[:, [0, 6]], test[:, [0, 6]])
    assert_allclose(actual[:, 4].astype(float), values, rtol=1e-13, atol=2e-15)
    for identity in np.unique(ids):
        rows = actual[ids == identity]
        assert len(set(rows[:, 0])) == 1
        positions = rows[:, 1].astype(float)
        first, last = target["centromeres"][rows[0, 0]]
        assert not ((positions < first).any() and (positions > last).any())


def test_arm_evidence_checksums():
    manifest = json.loads((EDGES / "manifest.json").read_text())
    for name, record in manifest["files"].items():
        assert hashlib.sha256((EDGES / name).read_bytes()).hexdigest() == record["sha256"]


def test_empty_declared_chromosome_fails_explicitly():
    target = load_manifest(FIXTURE / "inputs/target", "target")
    target["chromosomes"].append("chr3")
    target["centromeres"]["chr3"] = [7000, 11000]
    prepared = load_manifest(FIXTURE / "inputs/prepared", "prepared")
    with np.load(FIXTURE / "inputs/prepared" / prepared["samples"]["Test1"]) as data:
        matrix = data["matrix"]
    with pytest.raises(ValueError, match="empty or single-window chromosome: chr3"):
        segment_profile(matrix, np.linspace(-1, 1, len(matrix)), target, DEFAULTS)
