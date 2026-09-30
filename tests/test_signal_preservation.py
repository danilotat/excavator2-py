"""Original normalization and pooling contracts, without calibration inputs."""

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from excavator2.analyze import ratios
from excavator2.normalization import normalize

FIXTURE = Path(__file__).parent / "fixtures/legacy-analysis"


def matrix(values, classes):
    return np.column_stack(
        [
            np.repeat("chr1", len(values)),
            np.arange(len(values)),
            np.arange(len(values)),
            np.arange(len(values)) + 1,
            np.repeat("id", len(values)),
            values,
            classes,
        ]
    )


@pytest.mark.parametrize("experiment", ["paired", "pooling"])
def test_ratios_match_original_r_fixtures(experiment):
    import json

    folder = FIXTURE / "inputs/prepared"
    manifest = json.loads((folder / "manifest.json").read_text())
    samples = {}
    for name, filename in manifest["samples"].items():
        with np.load(folder / filename) as data:
            samples[name] = data["matrix"]
    for name in ["Test1", "Test2"]:
        controls = (
            [samples["Control" + name[-1]]]
            if experiment == "paired"
            else [samples["Control1"], samples["Control2"]]
        )
        expected = np.loadtxt(
            FIXTURE / "expected" / experiment / name / f"HSLMResults_{name}.txt",
            dtype=str,
            skiprows=1,
        )
        assert_allclose(
            ratios(samples[name], controls), expected[:, 4].astype(float), atol=2e-14, rtol=1e-13
        )


def test_pooling_uses_normalized_values_and_separate_class_medians():
    classes = ["IN"] * 3 + ["OUT"] * 3
    test = np.array([100, 150, 300, 200, 200, 400], dtype=float)
    controls = np.array([[100, 100, 100, 100, 100, 100], [200, 100, 300, 200, 200, 300]])
    expected = np.log2(test / controls.mean(axis=0))
    expected[:3] -= np.median(expected[:3])
    expected[3:] -= np.median(expected[3:])
    assert_allclose(ratios(matrix(test, classes), [matrix(c, classes) for c in controls]), expected)


def test_original_zero_replacement_uses_each_class_minimum():
    target = np.array(
        [
            ["chr1", str(i * 20), str(i * 20 + 10), f"a{i}", kind]
            for i, kind in enumerate(["IN"] * 3 + ["OUT"] * 3)
        ]
    )
    trace = normalize([0, 100, 200, 0, 300, 400], target, np.full(6, 0.5), np.full(6, 0.9))
    assert_allclose(trace["normalized"], [10, 10, 20, 30, 30, 40])
    assert_allclose(trace["weighted"], [0, 10, 20, 0, 30, 40])
