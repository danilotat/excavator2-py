"""Chromosome-keyed legacy feature imports and alignment validation."""

import json

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from excavator2.artifacts import load_preparation_exports


def export(path, values):
    path.mkdir(parents=True, exist_ok=True)
    values = np.asarray(values)
    character = values.dtype.kind in "US"
    metadata = {
        "order": "F",
        "endian": "little",
        "length": values.size,
        "dim": list(values.shape) if values.ndim > 1 else None,
        "type": "character" if character else "double",
    }
    (path / "metadata.json").write_text(json.dumps(metadata))
    np.zeros(values.size, dtype="<i4").tofile(path / "missing.bin")
    if character:
        (path / "values.json").write_text(json.dumps(values.flatten(order="F").tolist()))
    else:
        values.astype("<f8").flatten(order="F").tofile(path / "data.bin")


@pytest.fixture
def exported(tmp_path):
    export(
        tmp_path / "panel.RData/MyTarget",
        [
            ["chr1", "10", "20", "a1", "IN"],
            ["chr2", "30", "40", "a2", "OUT"],
            ["chr2", "50", "60", "a3", "IN"],
        ],
    )
    export(tmp_path / "GCC/GCC.chr1.RData/GCContent", [0.1])
    export(tmp_path / "GCC/GCC.chr2.RData/GCContent", [0.2, 0.3])
    export(tmp_path / "MAP/Map.chr1.RData/MapMed", [0.4])
    export(tmp_path / "MAP/Map.chr2.RData/MapMed", [0.5, 0.6])
    # Only one directory has an extra chromosome sorting before the requested ones.
    export(tmp_path / "GCC/GCC.chr0.RData/GCContent", [0.9])
    return tmp_path


def test_features_are_joined_by_chromosome(exported):
    actual = load_preparation_exports(exported, ["chr1", "chr2"])
    assert_array_equal(actual["gc"], [0.1, 0.2, 0.3])
    assert_array_equal(actual["mappability"], [0.4, 0.5, 0.6])


@pytest.mark.parametrize(
    "directory,prefix,obj", [("GCC", "GCC", "GCContent"), ("MAP", "Map", "MapMed")]
)
def test_missing_feature_does_not_select_another_chromosome(exported, directory, prefix, obj):
    path = exported / directory / f"{prefix}.chr2.RData"
    path.rename(path.with_name(f"{prefix}.chr3.RData"))
    with pytest.raises(ValueError, match=f"missing {directory} export for chr2"):
        load_preparation_exports(exported, ["chr1", "chr2"])


@pytest.mark.parametrize("values", [[0.1], [0.1, float("nan")], [[0.1, 0.2]]])
def test_invalid_features_fail_alignment_validation(exported, values):
    export(exported / "GCC/GCC.chr2.RData/GCContent", values)
    with pytest.raises(ValueError, match="GCC export does not align with chr2"):
        load_preparation_exports(exported, ["chr1", "chr2"])


def test_export_chromosome_order_must_match_target(exported):
    with pytest.raises(ValueError, match="declared chromosome order"):
        load_preparation_exports(exported, ["chr2", "chr1"])
