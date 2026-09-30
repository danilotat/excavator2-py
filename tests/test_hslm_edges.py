"""Retained legacy boundary evidence and corrected singleton behavior."""

import hashlib
import json
from pathlib import Path

import pytest
from numpy.testing import assert_array_equal

from excavator2.hslm import estimate_parameters, segment

FIXTURES = Path(__file__).parent / "fixtures/legacy-hslm/edges"


@pytest.mark.parametrize("backend", ["python", "native"])
def test_single_window_is_preserved(backend):
    parameters = estimate_parameters([0.1])
    result = segment([0.1], [100], parameters, backend=backend)
    assert_array_equal(result.filtered_breaks, [0, 1])
    assert_array_equal(result.values, [0.1])


def test_boundary_evidence_checksums():
    manifest = json.loads((FIXTURES / "manifest.json").read_text())
    for name, record in manifest["files"].items():
        assert hashlib.sha256((FIXTURES / name).read_bytes()).hexdigest() == record["sha256"]
