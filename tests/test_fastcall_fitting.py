"""Truncated objective and numerical-underflow regressions."""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy.integrate import quad
from scipy.special import logsumexp
from scipy.stats import norm, truncnorm

from excavator2 import _core
from excavator2.fastcall import fit_fastcall
from excavator2.reference import fastcall as reference


@pytest.fixture(params=["python", "native"])
def kernels(request):
    if request.param == "python":
        return reference.expectation, reference.maximization, reference.posterior
    return _core.fastcall_expectation, _core.fastcall_maximization, _core.fastcall_posterior


def parameters():
    return (
        np.array([-3.0, -1.0, 0.0, 0.58, 1.0]),
        np.full(5, 0.4),
        np.array([[-20.0, -1.3], [-1.3, -0.5], [-0.5, 0.35], [0.35, 0.9], [0.9, 20.0]]),
    )


def component_q(values, weights, sd, lower, upper):
    if np.isposinf(sd):
        return -weights.sum() * np.log(upper - lower)
    return np.dot(weights, truncnorm.logpdf(values, lower / sd, upper / sd, scale=sd))


def test_a20_counterexample_increases_truncated_objective(kernels):
    _, maximize, _ = kernels
    means, sd, bounds = parameters()
    values = np.array([0.2, 0.3])
    responsibilities = np.zeros((2, 5))
    responsibilities[:, 2] = 1
    fitted, priors = maximize(values, responsibilities, means, sd, bounds)
    before = component_q(values, np.ones(2), 0.4, -0.5, 0.35)
    after = component_q(values, np.ones(2), fitted[2], -0.5, 0.35)
    assert_allclose(fitted[2], 1.59228318, rtol=1e-7)
    assert_allclose(before, 0.2916493283198236, atol=1e-13)
    assert_allclose(after, 0.325201823605964, atol=1e-13)
    assert after > before
    assert_array_equal(priors, [0, 0, 1, 0, 0])
    assert_array_equal(fitted[[0, 1, 3, 4]], sd[[0, 1, 3, 4]])


@pytest.mark.parametrize(
    "lower, upper, values",
    [
        (-0.5, 0.35, [0.01, 0.02, 0.04]),
        (-0.5, 0.35, [0.2, 0.25]),
        (9.0, 10.0, [9.001, 9.01, 9.02]),
        (-10.0, -9.0, [-9.001, -9.01, -9.02]),
    ],
)
def test_mstep_stationarity_and_independent_second_moment(kernels, lower, upper, values):
    _, maximize, _ = kernels
    means, sd, bounds = parameters()
    bounds[2] = [lower, upper]
    values = np.array(values)
    weights = np.linspace(0.2, 1, len(values))
    responsibilities = np.zeros((len(values), 5))
    responsibilities[:, 2] = weights
    fitted, _ = maximize(values, responsibilities, means, sd, bounds)
    sigma = fitted[2]
    assert np.isfinite(sigma)

    def q(scale):
        return component_q(values, weights, scale, lower, upper)

    assert q(sigma) >= q(sd[2]) - 1e-10
    derivative = (q(sigma * np.exp(1e-5)) - q(sigma * np.exp(-1e-5))) / 2e-5
    assert abs(derivative) < 2e-7
    moment, _ = quad(
        lambda x: x * x * truncnorm.pdf(x, lower / sigma, upper / sigma, scale=sigma),
        lower,
        upper,
        epsabs=1e-10,
        epsrel=1e-10,
    )
    assert_allclose(moment, np.dot(weights, values**2) / weights.sum(), rtol=1e-9, atol=1e-10)


@pytest.mark.parametrize("values, uniform", [([0.0, 0.0], False), ([-0.5, 0.35], True)])
def test_variance_floor_and_uniform_limit(kernels, values, uniform):
    expect, maximize, _ = kernels
    means, sd, bounds = parameters()
    values = np.array(values)
    r = np.zeros((2, 5))
    r[:, 2] = 1
    fitted, priors = maximize(values, r, means, sd, bounds)
    assert fitted[2] == (np.inf if uniform else 0.001)
    probabilities, likelihood = expect(values, means, fitted, priors, bounds)
    assert_array_equal(probabilities, r)
    assert_allclose(likelihood, component_q(values, np.ones(2), fitted[2], -0.5, 0.35), atol=1e-12)


def test_a21_respects_support_despite_tiny_densities(kernels):
    expect, _, _ = kernels
    means, _, bounds = parameters()
    sd, priors = np.full(5, 0.001), np.full(5, 0.2)
    probabilities, likelihood = expect(np.array([0.3]), means, sd, priors, bounds)
    assert_array_equal(probabilities, [[0, 0, 1, 0, 0]])
    assert_allclose(likelihood, np.log(0.2) + truncnorm.logpdf(0.3, -500, 350, scale=0.001))


def test_a21_uses_deviations_and_priors_instead_of_nearest_mean(kernels):
    _, _, posterior = kernels
    means, _, _ = parameters()
    sd = np.array([0.01, 0.01, 0.001, 0.001, 0.01])
    priors = np.array([0.01, 0.02, 0.9, 0.06, 0.01])
    values = np.array([0.29])
    logs = norm.logpdf(values[:, None], loc=means, scale=sd) + np.log(priors)
    expected = np.exp(logs - logsumexp(logs, axis=1, keepdims=True))
    actual = posterior(values, means, sd, priors)
    assert_allclose(actual, expected, atol=1e-14)
    assert np.argmax(actual[0]) == 4


@pytest.mark.parametrize("sd_value", [0.001, 0.03, 0.3, 1.0, 100.0, 1e8])
def test_boundary_competition_across_dynamic_range(kernels, sd_value):
    expect, _, _ = kernels
    means, _, bounds = parameters()
    sd, priors = np.full(5, sd_value), np.array([0.1, 0.2, 1e-200, 0.3, 0.4])
    values = np.array([-0.5, 0.35, 0.9])
    actual, likelihood = expect(values, means, sd, priors, bounds)
    logs = np.full((len(values), 5), -np.inf)
    for j, (lower, upper) in enumerate(bounds):
        inside = (values >= lower) & (values <= upper)
        if not inside.any():
            continue
        # Integrate in original units for broad curves to avoid tiny CDF differences.
        if sd_value > 100:
            integral, _ = quad(
                lambda x: np.exp(-0.5 * ((x - means[j]) / sd_value) ** 2), lower, upper
            )
            logs[inside, j] = (
                np.log(priors[j])
                - 0.5 * ((values[inside] - means[j]) / sd_value) ** 2
                - np.log(integral)
            )
        else:
            logs[inside, j] = np.log(priors[j]) + truncnorm.logpdf(
                values[inside],
                (lower - means[j]) / sd_value,
                (upper - means[j]) / sd_value,
                loc=means[j],
                scale=sd_value,
            )
    shifted = logs - logs.max(axis=1, keepdims=True)
    expected = np.exp(shifted) / np.exp(shifted).sum(axis=1, keepdims=True)
    assert_allclose(actual, expected, rtol=2e-11, atol=2e-13)
    assert_allclose(likelihood, logsumexp(logs, axis=1).sum(), rtol=1e-11, atol=1e-11)


def test_zero_priors_and_impossible_responsibilities_fail_explicitly(kernels):
    expect, maximize, _ = kernels
    means, sd, bounds = parameters()
    values = np.array([0.0])
    with pytest.raises(FloatingPointError, match="model support"):
        expect(values, means, sd, np.zeros(5), bounds)
    with pytest.raises(ValueError, match="mass"):
        maximize(values, np.zeros((1, 5)), means, sd, bounds)
    r = np.zeros((1, 5))
    r[0, 4] = 1
    with pytest.raises(ValueError, match="support"):
        maximize(values, r, means, sd, bounds)


@pytest.mark.parametrize("backend", ["python", "native"])
@pytest.mark.parametrize("copies", [1, 2, 3])
def test_convergence_tracks_likelihood_for_counts_divisible_by_five(backend, copies):
    values = np.tile([-1.3, -0.5, 0.35, 0.9, 0.0], copies)
    fit = fit_fastcall(values, backend=backend)
    assert fit.converged
    assert fit.iterations > 2
    assert np.all(np.diff(fit.trace[:, -1]) >= -1e-10)
    # A further M-step cannot materially change a converged fit.
    sd, priors = reference.maximization(
        values, fit.posterior, fit.means, fit.deviations, fit.bounds
    )
    assert_allclose(1 / sd**2, 1 / fit.deviations**2, rtol=2e-8, atol=1e-8)
    assert_allclose(priors, fit.priors, atol=1e-8, rtol=0)


@pytest.mark.parametrize("lower", [9.0, 40.0])
@pytest.mark.parametrize("sigma", [0.001, 0.1, 1.0])
def test_extreme_tail_matches_scaled_quadrature(kernels, lower, sigma):
    expect, _, _ = kernels
    means = np.linspace(-1, 1, 5)
    sd = np.full(5, sigma)
    priors = np.array([0.05, 0.1, 0.7, 0.1, 0.05])
    bounds = np.tile([lower, lower + 1], (5, 1))
    values = lower + np.array([0.0, 1.0, 3.0]) * sigma**2 / lower
    logs = np.empty((3, 5))
    for j, mean in enumerate(means):
        distance = lower - mean
        scale = sigma**2 / distance
        integral, _ = quad(
            lambda t: np.exp(-t - 0.5 * (sigma / distance * t) ** 2),
            0,
            min(100, 1 / scale),
            epsabs=1e-13,
            epsrel=1e-13,
        )
        offsets = values - lower
        logs[:, j] = (
            np.log(priors[j])
            - np.log(scale * integral)
            - 0.5 * (offsets / sigma) * ((offsets + 2 * distance) / sigma)
        )
    totals = logsumexp(logs, axis=1)
    actual, likelihood = expect(values, means, sd, priors, bounds)
    assert_allclose(actual, np.exp(logs - totals[:, None]), rtol=1e-11, atol=1e-13)
    assert_allclose(likelihood, totals.sum(), rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("backend", ["python", "native"])
def test_fit_is_invariant_to_segment_order(backend):
    values = np.array([-3.0, -1.3, -0.5, 0.35, 0.9, 0.0, 0.1, 0.4, 2.0])
    order = np.random.default_rng(551).permutation(len(values))
    first = fit_fastcall(values, backend=backend)
    shuffled = fit_fastcall(values[order], backend=backend)
    assert first.converged and shuffled.converged
    assert_allclose(first.posterior[order], shuffled.posterior, rtol=1e-8, atol=1e-10)
    assert_allclose(first.priors, shuffled.priors, rtol=1e-8, atol=1e-10)
    assert_allclose(first.trace[-1, -1], shuffled.trace[-1, -1], rtol=1e-12, atol=1e-12)
