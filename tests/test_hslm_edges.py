"""Boundary contracts captured from the original pure-R HSLM wrappers."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2.hslm import estimate_parameters, filter_breaks, reconstruct, segment

FIXTURES = Path(__file__).parent / "fixtures/legacy-hslm/edges"


def records(filename):
    with (FIXTURES / filename).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def vector(value):
    return np.array([float(x) for x in value.split(",")])


@pytest.mark.parametrize("case", records("filter.tsv"), ids=lambda case: case["case"])
def test_filter_boundaries_match_original(case):
    filtered = filter_breaks(vector(case["breaks"]), int(case["min_windows"]))
    assert_array_equal(filtered, vector(case["filtered"]))
    assert_array_equal(reconstruct(vector(case["values"]), filtered), vector(case["reconstructed"]))


@pytest.mark.parametrize("case", records("parameters.tsv"), ids=lambda case: case["case"])
def test_parameter_degeneracy_matches_original_domain(case):
    values = vector(case["values"])
    if case["smu"] == "NA":
        # Original quantile trimming leaves fewer than two values: var returns NA.
        with pytest.raises(ValueError, match="two retained values"):
            estimate_parameters(values)
    elif float(case["smu"]) == 0:
        # Original passes zero deviations onwards; the port rejects this domain.
        with pytest.raises(ValueError, match="positive finite global variance"):
            estimate_parameters(values)
    else:
        result = estimate_parameters(values)
        assert_allclose(
            [result.mi, result.smu, result.sepsilon],
            [float(case[key]) for key in ("mi", "smu", "sepsilon")],
            rtol=2e-14,
            atol=0,
        )


@pytest.mark.parametrize("backend", ["python", "native"])
def test_single_window_wrapper_rejected_with_valid_global_parameters(backend):
    assert (FIXTURES / "single-window-error.txt").read_text().strip() == (
        "missing value where TRUE/FALSE needed"
    )
    parameters = estimate_parameters([-1, -0.5, 0, 0.5, 1])
    with pytest.raises(ValueError, match="single-window arm"):
        segment([0.1], [100], parameters, backend=backend)


def test_boundary_evidence_checksums():
    manifest = json.loads((FIXTURES / "manifest.json").read_text())
    for name, record in manifest["files"].items():
        assert hashlib.sha256((FIXTURES / name).read_bytes()).hexdigest() == record["sha256"]
