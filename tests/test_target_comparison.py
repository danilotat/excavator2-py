"""The target concordance gate must reject changed science, not just broken hashes."""

import importlib.util
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from excavator2.artifacts import digest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "compare_target", ROOT / "tools/oracle/compare_target.py"
)
comparator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparator)


@pytest.fixture
def targets(tmp_path):
    fixture = ROOT / "tests/fixtures/legacy-preparation/pipeline/target"
    old, new = tmp_path / "old", tmp_path / "new"
    shutil.copytree(fixture, old)
    shutil.copytree(fixture, new)
    return old, new


def test_target_comparison_accepts_identical_arrays(targets):
    assert comparator.compare(*targets)["passed"]


@pytest.mark.parametrize("field", ["gc", "mappability", "target", "reference", "centromeres"])
def test_target_comparison_rejects_changed_content_with_valid_checksums(targets, field):
    old, new = targets
    path = new / "manifest.json"
    manifest = json.loads(path.read_text())
    if field == "centromeres":
        manifest["centromeres"][manifest["chromosomes"][0]][0] += 1
    else:
        filename = (
            manifest["preparation"]
            if field != "reference"
            else next(iter(manifest["references"].values()))
        )
        with np.load(new / filename) as archive:
            arrays = {k: archive[k] for k in archive.files}
        if field == "target":
            arrays[field][[0, 1]] = arrays[field][[1, 0]]
        elif field == "reference":
            arrays["matrix"][0, 1] = "N"
        else:
            arrays[field][0] += 1e-6
        np.savez_compressed(new / filename, **arrays)
        manifest["files"][filename] = digest(new / filename)
    path.write_text(json.dumps(manifest))
    with pytest.raises((ValueError, AssertionError)):
        comparator.compare(old, new)
