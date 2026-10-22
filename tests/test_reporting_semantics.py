"""Evidence for A18/A25 reporting and P15/H07 model assessments."""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy.special import ndtr

from excavator2.fastcall import (
    assign_labels,
    call_probabilities,
    correct_cellularity,
    fit_fastcall,
    segment_standard_errors,
    window_noise,
)
from excavator2.writers import write_results


@pytest.mark.parametrize("chromosome", ["chr1", "chrX"])
@pytest.mark.parametrize("purity", [1, 0.5, 0.05])
def test_corrected_reporting_and_floor(tmp_path, purity, chromosome):
    # Includes deletion, gain and a mixture incompatible with nonnegative tumor CN.
    mixture = np.log2([1 - purity / 2, 1 + purity / 2, 0.001])
    corrected = correct_cellularity(mixture, purity)
    calls = assign_labels(fit_fastcall(corrected).posterior)
    probabilities = call_probabilities(mixture, np.full(3, 0.1), cellularity=purity)
    rows = np.array(
        [
            [chromosome, str(i), str(i), str(i + 1), str(x), str(x), "IN"]
            for i, x in zip([10, 20, 30], mixture, strict=True)
        ]
    )
    np.savez(tmp_path / "reference.npz", matrix=np.array([["10", "A"], ["20", "C"], ["30", "G"]]))
    target = {
        "assembly": "synthetic",
        "references": {chromosome: "reference.npz"},
        "files": {"reference.npz": "test"},
    }
    write_results(
        tmp_path,
        "T",
        rows,
        np.arange(3),
        np.arange(1, 4),
        mixture,
        corrected,
        calls,
        probabilities,
        target,
        tmp_path,
    )
    table = np.genfromtxt(
        tmp_path / "FastCallResults_T.txt", names=True, dtype=None, encoding="utf-8"
    )
    assert_allclose(table["MixtureSegment"], mixture)
    assert_allclose(table["CorrectedSegment"], corrected)
    assert_allclose(table["DECNF"][:2], [1, 3])
    assert_allclose(table["DECNF"][2], 0.002 if purity == 1 else 0.0625)
    assert_array_equal(table["DECN"], np.rint(table["DECNF"]))
    called = probabilities[np.arange(3), calls.labels + 2]
    assert_allclose(table["ProbCall"], called)
    assert_allclose(table["ProbAltered"], 1 - probabilities[:, 2])
    for kind in ["Window", "Region"]:
        text = (tmp_path / f"EXCAVATOR{kind}Call_T.vcf").read_text()
        assert "not absolute CN" in text
        assert "not calibrated event confidence" in text
        assert "ID=CN," not in text and "ID=CNF," not in text
        records = [line.split("\t") for line in text.splitlines() if not line.startswith("#")]
        values = [dict(zip(r[8].split(":"), r[9].split(":"), strict=True)) for r in records]
        assert_allclose([float(v["DECNF"]) for v in values], table["DECNF"])
        assert_allclose([float(v["DECN"]) for v in values], table["DECN"])
        assert all(v["GT"] == "." for v in values)
        assert_allclose([float(v["FCP"]) for v in values], np.round(called, 2))
        assert_allclose([float(v["FPA"]) for v in values], np.round(table["ProbAltered"], 2))


@pytest.mark.parametrize("backend", ["python", "native"])
def test_fit_posterior_ignores_support_and_noise(backend):
    # H07 reproducer: the diagnostic fit sees one median per segment.
    medians = np.array([-2, -1, 0, 0.6, 1.2, 0.6])
    fit = fit_fastcall(medians, backend=backend)
    assert_array_equal(fit.posterior, np.eye(5)[[0, 1, 2, 3, 4, 3]])
    # More segments change priors, but interior class membership stays certain.
    fragmented = fit_fastcall(np.r_[medians, np.full(20, 0.6)], backend=backend)
    assert fragmented.priors[3] > fit.priors[3]
    assert_array_equal(fragmented.posterior[:6], fit.posterior)


def test_call_probabilities_depend_on_support_noise_and_threshold_distance():
    noise = {"IN": 0.34, "OUT": 0.57}
    # One 4-window and one 200-window IN segment, both with median 0.5 (a gain).
    classes = np.array(["IN"] * 204)
    starts, ends = np.array([0, 4]), np.array([4, 204])
    errors = segment_standard_errors(classes, starts, ends, noise)
    assert_allclose(errors, np.sqrt(np.pi / 2) * 0.34 / np.sqrt([4, 200]))
    gain = call_probabilities([0.5, 0.5], errors)[:, 3]
    expected = ndtr((0.9 - 0.5) / errors) - ndtr((0.35 - 0.5) / errors)
    assert_allclose(gain, expected)
    assert gain[0] < 0.75 < 0.99 < gain[1]
    # Noisier OUT windows widen a mixed segment, and rows are distributions.
    mixed = segment_standard_errors(np.array(["IN", "OUT"] * 10), [0], [20], noise)
    assert mixed[0] > segment_standard_errors(np.array(["IN"] * 20), [0], [20], noise)[0]
    probabilities = call_probabilities([-2, -0.6, 0, 0.36, 3], np.full(5, 0.2))
    assert_allclose(probabilities.sum(axis=1), 1)
    assert probabilities[3, 2] > 0.4  # just past the gain threshold


@pytest.mark.parametrize("purity", [1, 0.5, 0.05])
def test_call_probabilities_split_at_cellularity_mapped_edges(purity):
    edges = np.array([-1.3, -0.5, 0.35, 0.9])
    mixture = np.log2(purity * np.exp2(edges) + 1 - purity)
    probabilities = call_probabilities(mixture, np.full(4, 1e-5), cellularity=purity)
    for i, row in enumerate(probabilities):
        assert_allclose(row[[i, i + 1]], 0.5, atol=1e-12)
    # A mixture below the correction floor is a deep loss with certainty.
    floor = call_probabilities([np.log2(0.001)], [0.01], cellularity=purity)
    assert_allclose(floor[0, 0], 1)


def test_window_noise_excludes_steps_and_separates_classes():
    rng = np.random.default_rng(7)
    classes = np.where(np.arange(4000) % 4 == 0, "OUT", "IN")
    sigma = np.where(classes == "IN", 0.3, 0.6)
    levels = np.repeat([0.0, 1.0, -1.0, 0.585], 1000)
    values = levels + sigma * rng.standard_normal(4000)
    starts, ends = np.arange(0, 4000, 1000), np.arange(1000, 4001, 1000)
    noise = window_noise(values, classes, starts, ends)
    assert_allclose([noise["IN"], noise["OUT"]], [0.3, 0.6], rtol=0.08)
    # Probabilities of a segment depend on its own windows, not on fragmentation elsewhere.
    errors = segment_standard_errors(classes, starts, ends, noise)
    split = np.r_[0, 500, starts[1:]], np.r_[500, ends]
    split_errors = segment_standard_errors(classes, *split, noise)
    assert_allclose(split_errors[2:], errors[1:])


def test_reference_bases_match_scientific_positions(tmp_path):
    # Row text renders round positions in R style ("1e+05"); FRB stores plain integers.
    mixture = np.log2([0.5, 1.5])
    calls = assign_labels(fit_fastcall(mixture).posterior)
    rows = np.array(
        [
            ["chr1", "1", "1e+05", "100100", str(mixture[0]), str(mixture[0]), "IN"],
            ["chr1", "2", "2e+05", "200100", str(mixture[1]), str(mixture[1]), "IN"],
        ]
    )
    np.savez(tmp_path / "reference.npz", matrix=np.array([["100000", "A"], ["200000", "C"]]))
    target = {
        "assembly": "synthetic",
        "references": {"chr1": "reference.npz"},
        "files": {"reference.npz": "test"},
    }
    write_results(
        tmp_path,
        "T",
        rows,
        np.arange(2),
        np.arange(1, 3),
        mixture,
        mixture,
        calls,
        call_probabilities(mixture, np.full(2, 0.1)),
        target,
        tmp_path,
    )
    text = (tmp_path / "EXCAVATORWindowCall_T.vcf").read_text()
    records = [line.split("\t") for line in text.splitlines() if not line.startswith("#")]
    assert [(r[1], r[3]) for r in records] == [("100000", "A"), ("200000", "C")]
