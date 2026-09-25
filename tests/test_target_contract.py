"""Original target geometry fixtures and exact Python-port comparisons."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal

FIXTURES = Path(__file__).parent / "fixtures/legacy-target"
MANIFEST = json.loads((FIXTURES / "manifest.json").read_text())


def target(name):
    with np.load(FIXTURES / name / "expected.npz", allow_pickle=False) as data:
        return data["target"]


@pytest.mark.parametrize("name", MANIFEST["cases"])
def test_captured_target_case_integrity(name):
    case = MANIFEST["cases"][name]
    for filename, record in case["files"].items():
        content = (FIXTURES / name / filename).read_bytes()
        assert len(content) == record["size"]
        assert hashlib.sha256(content).hexdigest() == record["sha256"]
    assert (case["status"] == 0) == case["expected_success"]
    if case["expected_success"]:
        rows = target(name)
        assert rows.shape == (case["rows"], 5)
        assert_array_equal(rows[:, 3], [f"a{i}" for i in range(1, len(rows) + 1)])
        prefix = "" if name == "bare_names" else "chr"
        chromosomes = [f"{prefix}{c}" for c in [*range(1, 23), "X"]]
        assert list(dict.fromkeys(rows[:, 0])) == chromosomes
        with np.load(FIXTURES / name / "expected.npz") as data:
            assert_array_equal(data["chromosomes"], chromosomes)
    else:
        assert not (FIXTURES / name / "expected.npz").exists()
        error = (FIXTURES / name / "stderr.log").read_text()
        expected = (
            "'to' must be a finite number"
            if name == "no_eligible_gap"
            else "subscript out of bounds"
        )
        assert expected in error
        assert "Execution halted" in error


def test_chromosome_naming_contract():
    standard = target("standard")
    bare = target("bare_names")
    assert_array_equal([f"chr{c}" for c in bare[:, 0]], standard[:, 0])
    assert_array_equal(bare[:, 1:], standard[:, 1:])
    assert_array_equal(target("coordinate_names_ignored"), standard)


def test_gap_filter_uses_endpoints_not_overlap():
    rows = target("gap_endpoints")
    inside = rows[(rows[:, 0] == "chr1") & (rows[:, 4] == "IN")]
    assert_array_equal(
        inside[:, 1:3], [["500", "560"], ["1200", "1700"], ["1600", "1650"], ["2000", "2100"]]
    )


def test_unsorted_bed_preserves_duplicate_start_order():
    rows = target("unsorted_overlap")
    inside = rows[(rows[:, 0] == "chr1") & (rows[:, 4] == "IN")]
    assert_array_equal(
        inside[:, 1:3], [["500", "560"], ["550", "700"], ["550", "650"], ["1300", "1400"]]
    )
    assert len(rows) != len(target("standard"))


def test_off_target_endpoint_and_flank_quirks():
    rows = target("standard")
    chr1 = rows[rows[:, 0] == "chr1"]
    chr2 = rows[rows[:, 0] == "chr2"]
    assert_array_equal(chr1[0, 1:3], ["201", "300"])
    assert_array_equal(chr2[0, 1:3], ["0", "99"])
    assert_array_equal(chr1[-1, 1:3], ["3101", "3200"])
    assert any(list(row[[1, 2, 4]]) == ["501", "600", "OUT"] for row in chr1)


@pytest.mark.parametrize("name", MANIFEST["cases"])
def test_python_geometry_matches_original(name):
    from excavator2.target import target_geometry

    folder = FIXTURES / name
    arguments = (folder / "target.bed", folder / "chromosomes.tsv", folder / "gaps.tsv", 100)
    if MANIFEST["cases"][name]["expected_success"]:
        actual = target_geometry(*arguments)
        assert_array_equal(actual, target(name))
        with np.load(folder / "expected.npz") as data:
            assert_array_equal(list(dict.fromkeys(actual[:, 0])), data["chromosomes"])
    else:
        errors = {
            "no_eligible_gap": "no eligible off-target gap for chr1",
            "no_alternate_gap": "requires an alternate-contig gap row",
            "missing_chromosome_gap": "requires gaps for chr2",
        }
        with pytest.raises(ValueError, match=errors[name]):
            target_geometry(*arguments)


@pytest.mark.parametrize("window", [0, -1, 9, 10.5, True])
def test_geometry_rejects_invalid_window(window):
    from excavator2.target import target_geometry

    folder = FIXTURES / "standard"
    with pytest.raises(ValueError, match="window must be an integer of at least 10"):
        target_geometry(
            folder / "target.bed", folder / "chromosomes.tsv", folder / "gaps.tsv", window
        )


def test_geometry_ignores_extra_bed_columns(tmp_path):
    from excavator2.target import target_geometry

    folder = FIXTURES / "standard"
    bed = tmp_path / "annotated.bed"
    bed.write_text(
        "".join(
            f"{line}\tregion{i}\n"
            for i, line in enumerate((folder / "target.bed").read_text().splitlines())
        )
    )
    actual = target_geometry(bed, folder / "chromosomes.tsv", folder / "gaps.tsv", 100)
    assert_array_equal(actual, target("standard"))


def test_geometry_rejects_incomplete_coordinate_table(tmp_path):
    from excavator2.target import target_geometry

    folder = FIXTURES / "standard"
    coordinates = tmp_path / "chromosomes.tsv"
    coordinates.write_text("chr1\t0\t3000\n")
    with pytest.raises(ValueError, match="23 canonical coordinate rows"):
        target_geometry(folder / "target.bed", coordinates, folder / "gaps.tsv", 100)
