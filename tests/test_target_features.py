"""Exact reference-feature comparisons against unchanged legacy shell/R tools."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from excavator2.target_features import target_features

FIXTURES = Path(__file__).parent / "fixtures/legacy-target-features"


def test_feature_fixture_integrity():
    manifest = json.loads((FIXTURES / "manifest.json").read_text())
    for name, record in manifest["files"].items():
        content = (FIXTURES / name).read_bytes()
        assert len(content) == record["size"]
        assert hashlib.sha256(content).hexdigest() == record["sha256"]


@pytest.mark.parametrize("name", ["standard", "bare_names", "past_end"])
def test_features_match_original(name):
    target = np.loadtxt(FIXTURES / name / "target.tsv", dtype=str, ndmin=2)
    actual = target_features(target, FIXTURES / "reference.fa", FIXTURES / "reference.bw")
    with np.load(FIXTURES / name / "expected.npz", allow_pickle=False) as expected:
        assert set(expected.files) == {
            f"{kind}_{chrom}" for chrom in actual for kind in ["GCC", "MAP", "FRB"]
        }
        for chrom, features in actual.items():
            for kind, key in [("GCC", "gc"), ("MAP", "mappability"), ("FRB", "first_base")]:
                assert_array_equal(features[key], expected[f"{kind}_{chrom}"])


@pytest.mark.parametrize("name", ["zero_start", "last_base", "missing_contig", "partial_skips"])
def test_feature_failures_match_original(name):
    manifest = json.loads((FIXTURES / "manifest.json").read_text())
    assert manifest["cases"][name]["status"] != 0
    target = np.loadtxt(FIXTURES / name / "target.tsv", dtype=str, ndmin=2)
    message = (
        "positive valid coordinates"
        if name == "zero_start"
        else "first-base extraction has no rows"
    )
    with pytest.raises(ValueError, match=message):
        target_features(target, FIXTURES / "reference.fa", FIXTURES / "reference.bw")
