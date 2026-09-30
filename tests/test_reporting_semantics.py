"""Evidence for A18/A25 reporting and P15/H07 model assessments."""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2.fastcall import assign_labels, correct_cellularity, fit_fastcall
from excavator2.writers import write_results


@pytest.mark.parametrize("chromosome", ["chr1", "chrX"])
@pytest.mark.parametrize("purity", [1, 0.5, 0.05])
def test_corrected_reporting_and_floor(tmp_path, purity, chromosome):
    # Includes deletion, gain and a mixture incompatible with nonnegative tumor CN.
    mixture = np.log2([1 - purity / 2, 1 + purity / 2, 0.001])
    corrected = correct_cellularity(mixture, purity)
    calls = assign_labels(fit_fastcall(corrected).posterior)
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


@pytest.mark.parametrize("backend", ["python", "native"])
def test_probability_is_class_membership_not_support_confidence(backend):
    medians = np.array([-2, -1, 0, 0.6, 1.2, 0.6])
    fit = fit_fastcall(medians, backend=backend)
    assert_array_equal(fit.posterior, np.eye(5)[[0, 1, 2, 3, 4, 3]])
    # More segments change priors, but interior class membership stays certain.
    fragmented = fit_fastcall(np.r_[medians, np.full(20, 0.6)], backend=backend)
    assert fragmented.priors[3] > fit.priors[3]
    assert_array_equal(fragmented.posterior[:6], fit.posterior)
