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


def test_geometry_sorts_merges_and_bounds_windows(tmp_path):
    from excavator2.target import target_geometry

    folder = FIXTURES / "unsorted_overlap"
    args = (folder / "chromosomes.tsv", folder / "gaps.tsv", 100)
    actual = target_geometry(folder / "target.bed", *args)
    bed = tmp_path / "reversed.bed"
    bed.write_text("\n".join(reversed((folder / "target.bed").read_text().splitlines())))
    assert_array_equal(actual, target_geometry(bed, *args))
    inside = actual[(actual[:, 0] == "chr1") & (actual[:, 4] == "IN")]
    assert_array_equal(inside[:, 1:3], [["500", "700"], ["1300", "1400"]])
    for chrom in set(actual[:, 0]):
        rows = actual[actual[:, 0] == chrom]
        positions = rows[:, 1:3].astype(int)
        assert np.all(positions[:, 1] <= 3000)
        for left, right in positions[rows[:, 4] == "OUT"]:
            assert right - left + 1 == 100
            for a, b in positions[rows[:, 4] == "IN"]:
                assert right <= a - 200 or left > b + 200


def test_geometry_bounds_follow_names(tmp_path):
    from excavator2.target import target_geometry

    folder = FIXTURES / "standard"
    coordinates = tmp_path / "coordinates.tsv"
    lines = (folder / "chromosomes.tsv").read_text().splitlines()
    lines[0] = "chr1\t0\t2000"
    coordinates.write_text("\n".join(reversed(lines)))
    bed = tmp_path / "canonical.bed"
    bed.write_text("chr1\t500\t560\nchr1\t1300\t1400\n")
    actual = target_geometry(bed, coordinates, folder / "gaps.tsv", 100)
    assert actual[actual[:, 0] == "chr1", 2].astype(int).max() <= 2000
    assert actual[actual[:, 0] == "chrX", 2].astype(int).max() == 3000


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
            if not line.startswith("chrY\t")
        )
    )
    plain = tmp_path / "plain.bed"
    plain.write_text(
        "\n".join("\t".join(row.split("\t")[:3]) for row in bed.read_text().splitlines())
    )
    actual = target_geometry(bed, folder / "chromosomes.tsv", folder / "gaps.tsv", 100)
    assert_array_equal(
        actual,
        target_geometry(plain, folder / "chromosomes.tsv", folder / "gaps.tsv", 100),
    )


def test_geometry_rejects_incomplete_coordinate_table(tmp_path):
    from excavator2.target import target_geometry

    folder = FIXTURES / "standard"
    coordinates = tmp_path / "chromosomes.tsv"
    coordinates.write_text("chr1\t0\t3000\n")
    with pytest.raises(ValueError, match="23 canonical coordinate rows"):
        target_geometry(folder / "target.bed", coordinates, folder / "gaps.tsv", 100)


@pytest.mark.parametrize("kind", ["duplicate", "missing"])
def test_geometry_rejects_ambiguous_named_bounds(tmp_path, kind):
    from excavator2.target import target_geometry

    folder = FIXTURES / "standard"
    rows = (folder / "chromosomes.tsv").read_text().splitlines()
    rows[-1] = "chr1\t0\t3000" if kind == "duplicate" else "chrY\t0\t3000"
    coords = tmp_path / "coordinates.tsv"
    coords.write_text("\n".join(rows))
    with pytest.raises(ValueError, match="duplicate|all 23 canonical"):
        target_geometry(folder / "target.bed", coords, folder / "gaps.tsv", 100)


def test_target_without_eligible_out_gap_keeps_inside(tmp_path):
    from excavator2.target import target_geometry

    folder = FIXTURES / "no_eligible_gap"
    gaps = tmp_path / "gaps.tsv"
    gaps.write_text((folder / "gaps.tsv").read_text().replace("chr1\t", "chr1_unused\t"))
    with gaps.open("a") as handle:
        handle.write("0\tchr1\t3000\t3001\n")
    actual = target_geometry(folder / "target.bed", folder / "chromosomes.tsv", gaps, 100)
    assert set(actual[actual[:, 0] == "chr1", 4]) == {"IN"}


@pytest.mark.parametrize("contig", ["chrY", "chrM", "MT", "chrUn", "1"])
def test_unsupported_bed_contigs_fail_instead_of_disappearing(tmp_path, contig):
    from excavator2.target import target_geometry

    folder = FIXTURES / "standard"
    bed = tmp_path / "target.bed"
    bed.write_text(f"chr1\t500\t560\n{contig}\t800\t850\n")
    with pytest.raises(ValueError, match=f"unsupported BED contigs: {contig}"):
        target_geometry(bed, folder / "chromosomes.tsv", folder / "gaps.tsv", 100)


@pytest.mark.parametrize(
    "interval,kept",
    [
        ((400, 700), False),
        ((520, 550), False),
        ((450, 500), False),
        ((600, 650), True),
        ((400, 499), True),
        ((550, 650), False),
    ],
)
def test_gap_filter_excludes_full_interval_overlap(tmp_path, interval, kept):
    from excavator2.target import target_geometry

    folder = FIXTURES / "standard"
    bed = tmp_path / "target.bed"
    bed.write_text(f"chr1\t{interval[0]}\t{interval[1]}\n")
    gaps = tmp_path / "gaps.tsv"
    gaps.write_text((folder / "gaps.tsv").read_text() + "0\tchr1\t500\t600\n")
    actual = target_geometry(bed, folder / "chromosomes.tsv", gaps, 100)
    inside = actual[(actual[:, 0] == "chr1") & (actual[:, 4] == "IN")]
    assert bool(len(inside)) == kept
