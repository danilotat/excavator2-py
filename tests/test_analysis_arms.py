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

REJECTED = {"no_short_inside", "endpoints", "inside"}


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["case"])
@pytest.mark.parametrize("sample", ["Test1", "Test2"])
def test_consistent_arm_geometry(case, sample):
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
    if case["case"] in REJECTED:
        with pytest.raises(ValueError, match="overlap centromere"):
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
    with pytest.raises(ValueError, match="empty chromosome: chr3"):
        segment_profile(matrix, np.linspace(-1, 1, len(matrix)), target, DEFAULTS)


@pytest.mark.parametrize("interval", [(90, 110), (90, 100), (200, 210), (50, 250)])
@pytest.mark.parametrize("has_short_arm", [False, True])
def test_any_centromere_overlap_is_rejected(interval, has_short_arm):
    windows = [(10, 20)] if has_short_arm else []
    windows += [interval, (300, 320)]
    matrix = np.array(
        [["chr1", str((a + b) / 2), str(a), str(b), "id", "1", "IN"] for a, b in windows]
    )
    target = {"chromosomes": ["chr1"], "centromeres": {"chr1": [100, 200]}}
    with pytest.raises(ValueError, match="overlap centromere"):
        segment_profile(matrix, np.zeros(len(matrix)), target, DEFAULTS)


def test_singleton_chromosome_has_a_segment():
    matrix = np.array([["chr1", "10", "5", "15", "id", "1", "IN"]])
    target = {"chromosomes": ["chr1"], "centromeres": {"chr1": [100, 200]}}
    rows, path, indices, ids = segment_profile(matrix, np.array([0.5]), target, DEFAULTS)
    assert_array_equal(indices, [0])
    assert_array_equal(ids, [0])
    assert float(rows[0, 5]) == 0.5
