"""Regression evidence for P14, P13 and A17 before segmentation/calling."""

import json

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2.analyze import calibration_policy, ratios, run_analysis
from excavator2.normalization import normalize


def profile(test, controls, test_exposure=1, control_exposures=None):
    controls = np.atleast_2d(controls)
    return ratios(
        test,
        controls,
        test_exposure=test_exposure,
        control_exposures=control_exposures or [1] * len(controls),
    )


def test_zero_depth_is_local_and_depth_aware():
    test = np.r_[np.full(90, 100), np.zeros(9), 100]
    control = np.full(100, 100)
    before = profile(test, control)
    test[-1] = 1
    after = profile(test, control)
    assert_array_equal(before[:-1], after[:-1])
    assert_allclose(before[90:99], np.log2(0.5 / 100.5))
    assert profile([0], [1000])[0] < profile([0], [100])[0] < 0
    # An all-zero test class remains a finite loss; no sample minimum is needed.
    assert np.isfinite(profile(np.zeros(100), control)).all()


@pytest.mark.parametrize("test,control", [([0], [0]), ([100], [0]), ([0, 0], [0, 0])])
def test_reference_zero_is_unsupported(test, control):
    with pytest.raises(ValueError, match="zero-depth reference"):
        profile(test, control)
    with pytest.raises(ValueError, match="zero-depth reference"):
        profile(test, [control, np.ones(len(test))])


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -1, 1.5])
def test_missing_or_invalid_depth_is_not_a_zero(invalid):
    with pytest.raises(ValueError, match="raw counts"):
        profile([invalid], [100])
    with pytest.raises(ValueError, match="raw counts"):
        profile([100], [invalid])


@pytest.mark.parametrize("covariate", ["gc", "map", "size"])
def test_covariate_correlated_gain_and_shared_technical_bias(covariate):
    n = 100
    altered = np.arange(n) >= 80
    length = np.where(altered, 60, 40) if covariate == "size" else np.full(n, 40)
    feature = np.where(altered, 0.62, 0.42)
    gc = feature if covariate == "gc" else np.full(n, 0.5)
    mappability = feature if covariate == "map" else np.full(n, 0.9)
    target = np.column_stack(
        [np.repeat("chr1", n), np.zeros(n, int), length, np.repeat("gene", n), np.repeat("IN", n)]
    )
    # Same technical response in both samples, and a 1.5x biological effect in one bin.
    control = length * np.where(altered, 200, 100)
    test = control * np.where(altered, 1.5, 1)
    trace = normalize(test, target, gc, mappability)
    assert_allclose(trace["normalized"], test / length)
    retained = profile(test, control)
    assert_allclose(retained[altered], np.log2(1.5), atol=0.0001)
    assert_allclose(retained[~altered], 0, atol=1e-14)
    assert_allclose(profile(control, control), 0, atol=1e-14)


@pytest.mark.parametrize("altered_class", ["IN", "OUT"])
def test_majority_altered_class_does_not_redefine_baseline(altered_class):
    classes = np.repeat(["IN", "OUT"], 100)
    altered = (classes == altered_class) & (np.arange(200) % 100 < 60)
    values = profile(np.where(altered, 150, 100), np.full(200, 100))
    assert_allclose(values[altered], np.log2(150.5 / 100.5))
    assert_allclose(values[~altered], 0, atol=1e-14)


def test_whole_profile_scaling_uses_external_exposure_not_median():
    assert_allclose(profile([200] * 100, [100] * 100), np.log2(200.5 / 100.5))
    # Doubling sampling exposure is technical scaling, not a genome-wide gain.
    assert_allclose(profile([200] * 100, [100] * 100, test_exposure=2), 0, atol=0.004)
    # Controls have equal weight only after their different exposures are removed.
    pooled = profile([150], [[100], [200]], control_exposures=[1, 2])
    assert_allclose(pooled, np.log2(150.5 / ((100.5 + 200.5 / 2) / 2)))


@pytest.mark.parametrize(
    "patch",
    [
        {"baseline": "unknown"},
        {"baseline": "haploid-X-reference"},
        {"baseline": "unknown-PAR-reference"},
        {"baseline": "triploid-reference"},
        {"bias": "sample-specific"},
        {"exposure_source": "profile-median"},
        {"exposures": {"T": 0, "C": 1}},
        {"exposures": {"T": True, "C": 1}},
        {"exposures": {"T": float("inf"), "C": 1}},
        {"exposures": {"T": 1}},
    ],
)
def test_unsupported_baseline_and_bias_are_rejected(tmp_path, patch):
    policy = {
        "baseline": "diploid-reference",
        "bias": "shared",
        "exposure_source": "independent",
        "exposures": {"T": 1, "C": 1},
    }
    policy.update(patch)
    path = tmp_path / "calibration.json"
    path.write_text(json.dumps(policy))
    with pytest.raises(ValueError):
        calibration_policy(path, ["T", "C"])


def test_no_implicit_baseline():
    with pytest.raises(ValueError, match="requires --calibration"):
        calibration_policy(None, ["T", "C"])


def test_raw_policy_without_counts_cannot_reach_calling(tmp_path, monkeypatch):
    import shutil
    from pathlib import Path

    from excavator2 import analyze

    fixture = Path(__file__).parent / "fixtures/legacy-analysis"
    prepared = tmp_path / "prepared"
    shutil.copytree(fixture / "inputs/prepared", prepared)
    manifest_path = prepared / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["depth_policy"] = "raw-counts-v1"
    manifest_path.write_text(json.dumps(manifest))
    calibration = tmp_path / "calibration.json"
    calibration.write_text(
        json.dumps(
            {
                "baseline": "diploid-reference",
                "bias": "shared",
                "exposure_source": "independent",
                "exposures": dict.fromkeys(manifest["samples"], 1),
            }
        )
    )
    monkeypatch.setattr(analyze, "fit_fastcall", lambda *a, **kw: pytest.fail("reached calling"))
    with pytest.raises(ValueError, match="lacks raw counts"):
        run_analysis(
            fixture / "samples.yaml",
            prepared,
            fixture / "inputs/target",
            tmp_path / "results",
            "paired",
            calibration=calibration,
        )
    assert not (tmp_path / "results").exists()
