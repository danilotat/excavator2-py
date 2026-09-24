"""Native boundary safety and independent kernel comparisons."""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2 import _core
from excavator2.fastcall import fit_fastcall
from excavator2.reference import fastcall as reference


def parameters():
    means = np.array([-3.0, -1.0, 0.0, 0.58, 1.0])
    deviations = np.array([0.15, 0.2, 0.08, 0.12, 0.3])
    priors = np.array([0.05, 0.1, 0.7, 0.1, 0.05])
    edges = np.array([-20.0, -1.3, -0.5, 0.35, 0.9, 20.0])
    bounds = np.column_stack([edges[:-1], edges[1:]])
    return means, deviations, priors, bounds


@pytest.mark.parametrize("size", [1, 17, 10003])
def test_kernels_match_python_and_leave_inputs_unchanged(size):
    rng = np.random.default_rng(731)
    values = np.ascontiguousarray(rng.normal(size=size))
    means, deviations, priors, bounds = parameters()
    inputs = [values, means, deviations, priors, bounds]
    original = [array.copy() for array in inputs]
    for array in inputs:
        array.flags.writeable = False
    actual = _core.fastcall_expectation(values, means, deviations, priors, bounds)
    expected = reference.expectation(values, means, deviations, priors, bounds)
    assert_allclose(actual, expected, rtol=1e-13, atol=2e-15)
    actual_sd, actual_prior = _core.fastcall_maximization(values, actual, means, deviations)
    expected_sd, expected_prior = reference.maximization(values, expected, means, deviations)
    assert_allclose(actual_sd, expected_sd, rtol=1e-13, atol=1e-14)
    assert_allclose(actual_prior, expected_prior, rtol=1e-13, atol=1e-14)
    assert_allclose(
        _core.fastcall_posterior(values, means, deviations, priors),
        reference.posterior(values, means, deviations, priors),
        rtol=1e-13,
        atol=2e-15,
    )
    for array, before in zip(inputs, original, strict=True):
        assert_array_equal(array, before)


def test_underflow_uses_first_nearest_component():
    means = np.array([-3.0, -1.0, 0.0, 0.58, 1.0])
    deviations = np.full(5, 0.001)
    priors = np.full(5, 0.2)
    values = np.array([-2.0, -100.0, 100.0])
    result = _core.fastcall_posterior(values, means, deviations, priors)
    assert_array_equal(result, np.eye(5)[[0, 0, 4]])


def test_truncated_infinity_and_nan_follow_legacy_fallback():
    # Both upper-tail CDFs round to one: inside density/0 -> Inf -> 100;
    # outside 0/0 remains NaN. Do not "fix" this with a stable tail formula.
    values = np.array([9.5, 0.0])
    means = np.zeros(5)
    deviations = np.ones(5)
    priors = np.full(5, 0.2)
    bounds = np.tile([9.0, 10.0], (5, 1))
    actual = _core.fastcall_expectation(values, means, deviations, priors, bounds)
    expected = reference.expectation(values, means, deviations, priors, bounds)
    assert_array_equal(actual, expected)
    assert_array_equal(actual[0], np.full(5, 0.2))
    assert actual[1, 0] == 1
    assert np.isnan(actual[1, 1:]).all()


@pytest.mark.parametrize("values", [np.arange(6.0)[::2], np.arange(3, dtype=np.float32), [0.0]])
def test_native_boundary_rejects_implicit_copies(values):
    means, deviations, priors, _ = parameters()
    with pytest.raises(TypeError):
        _core.fastcall_posterior(values, means, deviations, priors)


def test_native_boundary_rejects_unaligned_buffer():
    means, deviations, priors, _ = parameters()
    values = np.ndarray((2,), dtype=np.float64, buffer=bytearray(17), offset=1)
    with pytest.raises(ValueError, match="aligned"):
        _core.fastcall_posterior(values, means, deviations, priors)


@pytest.mark.parametrize(
    "case",
    [
        "empty",
        "rank",
        "means",
        "sd_zero",
        "sd_nan",
        "negative_prior",
        "bounds",
        "unordered",
        "weights",
        "negative_weights",
    ],
)
def test_native_boundary_rejects_invalid_inputs(case):
    values = np.array([-1.0, 0.0, 1.0])
    means, deviations, priors, bounds = parameters()
    weights = np.ones((3, 5))
    if case == "empty":
        values = np.array([])
    elif case == "rank":
        values = values.reshape(1, 3)
    elif case == "means":
        means = means[:-1].copy()
    elif case == "sd_zero":
        deviations[0] = 0
    elif case == "sd_nan":
        deviations[0] = np.nan
    elif case == "negative_prior":
        priors[0] = -1
    elif case == "bounds":
        bounds = np.zeros((5, 1))
    elif case == "unordered":
        bounds[0] = [1, -1]
    elif case == "weights":
        weights = np.ones((3, 4))
    else:
        weights[0, 0] = -1
    with pytest.raises(ValueError):
        if case in ("weights", "negative_weights"):
            _core.fastcall_maximization(values, weights, means, deviations)
        else:
            _core.fastcall_expectation(values, means, deviations, priors, bounds)


def test_high_level_api_accepts_strided_values_and_rejects_unknown_backend():
    values = np.linspace(-3, 1, 34)[::2]
    native = fit_fastcall(values)
    python = fit_fastcall(values, backend="python")
    assert native.iterations == python.iterations
    assert_allclose(native.posterior, python.posterior, rtol=1e-13, atol=2e-15)
    with pytest.raises(ValueError, match="backend"):
        fit_fastcall(values, backend="typo")


def test_native_calls_have_no_shared_mutable_state():
    means, deviations, priors, _ = parameters()
    values = np.linspace(-3, 1, 10003)

    def calculate(_):
        return _core.fastcall_posterior(values, means, deviations, priors)

    expected = calculate(0)
    with ThreadPoolExecutor(max_workers=2) as pool:
        for result in pool.map(calculate, range(4)):
            assert_array_equal(result, expected)
