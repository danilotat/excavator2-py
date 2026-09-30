"""Regression probes for H04, A22, A24, A23/A04 and A05."""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2.analyze import DEFAULTS, segment_profile, summarize
from excavator2.hslm import Parameters, estimate_parameters, filter_breaks, segment


def fit(values, *, parameters=None, positions=None, backend="native"):
    values = np.asarray(values)
    return segment(
        values,
        np.arange(len(values)) * 1000 if positions is None else positions,
        estimate_parameters(values) if parameters is None else parameters,
        min_windows=4,
        backend=backend,
    )


@pytest.mark.parametrize("sign", [-1, 1])
def test_extreme_plateaus_keep_boundary_and_original_medians(sign):
    values = sign * np.repeat([0.0, 1.2, 2.4, 0.0], [300, 80, 80, 300])
    values += np.tile([-0.01, 0.01], len(values) // 2)
    result = fit(values)
    assert_array_equal(result.filtered_breaks, [0, 300, 380, 460, 760])
    assert_allclose(result.values[[0, 300, 380, 460]], sign * np.array([0, 1.2, 2.4, 0]))


def test_high_precision_low_amplitude_steps():
    values = np.repeat([0, 0.31, 0.34, 0], [100, 80, 80, 100])
    result = fit(values, parameters=Parameters(0, 0.3, 0.002))
    assert_array_equal(result.filtered_breaks, [0, 100, 180, 260, 360])


def matrix_for(values, chromosomes=None):
    n = len(values)
    positions = np.arange(n) * 1000 + 1
    return np.array(
        [
            chromosomes if chromosomes is not None else ["chr1"] * n,
            positions,
            positions,
            positions + 100,
            values,
            values,
            ["IN"] * n,
        ]
    ).T


def test_distant_burden_does_not_change_local_event():
    local = np.repeat([0, 0.5, 0], [60, 15, 60]) + np.resize([-0.02, 0.02], 135)
    target = {"chromosomes": ["chr1", "chr2"], "centromeres": {"chr1": [-2, -1], "chr2": [-2, -1]}}
    results = []
    for burden in [0, 1, 2]:
        background = np.repeat([-burden, burden], 500) + np.tile([-0.02, 0.02], 500)
        values = np.r_[local, background]
        matrix = matrix_for(values, ["chr1"] * len(local) + ["chr2"] * 1000)
        rows, _, _, ids = segment_profile(matrix, values, target, DEFAULTS)
        results.append((rows[:135, 5], ids[:135]))
        starts, ends, _ = summarize(rows, ids)
        assert_array_equal(starts[:3], [0, 60, 75])
        assert_array_equal(ends[:3], [60, 75, 135])
    for result in results[1:]:
        assert_array_equal(result[0], results[0][0])
        assert_array_equal(result[1], results[0][1])


@pytest.mark.parametrize("amplitude", [0.1, 0.7])
def test_terminal_events_reverse_with_uneven_distances(amplitude):
    values = np.r_[np.full(5, amplitude), np.zeros(120)] + np.resize([-0.01, 0.01], 125)
    positions = np.cumsum(np.resize([100, 3000, 10000], len(values)))
    forward = fit(values, positions=positions)
    reverse = fit(values[::-1], positions=positions[-1] - positions[::-1])
    assert_allclose(forward.values, reverse.values[::-1])
    assert_array_equal(forward.filtered_breaks, len(values) - reverse.filtered_breaks[::-1])


@pytest.mark.parametrize("breaks", [[0, 3], [0, 3, 13], [0, 10, 13], [0, 10, 13, 23]])
def test_filter_preserves_endpoints_and_reverses(breaks):
    values = np.repeat(np.arange(len(breaks) - 1), np.diff(breaks)).astype(float)
    result = filter_breaks(breaks, 4, values)
    reverse = filter_breaks(len(values) - np.array(breaks)[::-1], 4, values[::-1])
    assert_array_equal(result, len(values) - reverse[::-1])
    assert result[0] == 0 and result[-1] == len(values)


def test_short_gain_merges_toward_similar_neighbor():
    values = np.repeat([-1, 1, 1.1], [10, 3, 10])
    assert_array_equal(filter_breaks([0, 10, 13, 23], 4, values), [0, 10, 23])


def test_support_counts_in_only_and_keeps_exact_threshold():
    values = np.repeat([-1, 1, 1.1], [10, 4, 10])
    breaks = [0, 10, 14, 24]
    support = np.ones(24, dtype=bool)
    assert_array_equal(filter_breaks(breaks, 4, values, support), breaks)
    support[11] = False
    assert_array_equal(filter_breaks(breaks, 4, values, support), [0, 10, 24])


def test_equal_medians_keep_arm_and_chromosome_identity():
    values = np.full(40, 0.58)
    matrix = matrix_for(values, ["chr1"] * 20 + ["chr2"] * 20)
    target = {
        "chromosomes": ["chr1", "chr2"],
        "centromeres": {"chr1": [9500, 9600], "chr2": [29500, 29600]},
    }
    rows, _, _, ids = segment_profile(matrix, values, target, DEFAULTS)
    starts, ends, medians = summarize(rows, ids)
    assert_array_equal(starts, [0, 10, 20, 30])
    assert_array_equal(ends, [10, 20, 30, 40])
    assert_allclose(medians, 0.58)


def test_class_noise_and_native_reference_agree():
    rng = np.random.default_rng(42)
    classes = np.resize(["IN", "OUT"], 120)
    values = rng.normal(size=120) * np.where(classes == "IN", 0.01, 0.1)
    p = estimate_parameters(values, classes=classes)
    assert p.deviations[1] > 5 * p.deviations[0]
    native = segment(values, np.arange(120), p, trace=True)
    python = segment(values, np.arange(120), p, trace=True, backend="python")
    assert_array_equal(native.path, python.path)
    assert_allclose(native.emissions, python.emissions, rtol=2e-13, atol=2e-13)
    assert_allclose(native.scores, python.scores, rtol=2e-13, atol=2e-13)


def test_neutral_noise_has_no_spurious_events():
    rng = np.random.default_rng(123)
    for sd in [0.01, 0.05, 0.2]:
        result = fit(rng.normal(0, sd, 1000))
        assert_array_equal(result.filtered_breaks, [0, 1000])


@pytest.mark.parametrize("case", ["short", "zero", "nan", "stride", "dtype", "alignment"])
def test_native_rejects_invalid_precision_buffers(case):
    from excavator2 import _core

    noise = np.ones(4)
    if case == "short":
        noise = noise[:3]
    elif case == "zero":
        noise[0] = 0
    elif case == "nan":
        noise[0] = np.nan
    elif case == "stride":
        noise = np.ones(8)[::2]
    elif case == "dtype":
        noise = noise.astype(np.float32)
    else:
        noise = np.ndarray((4,), dtype=np.float64, buffer=bytearray(33), offset=1)
        noise[:] = 1
    with pytest.raises((ValueError, TypeError)):
        _core.hslm_segment(
            np.zeros(4), np.array([0.0]), 0, 1, 1, np.full(3, 0.1), np.zeros(1), False, noise
        )
