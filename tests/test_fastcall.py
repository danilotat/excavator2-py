"""FastCall fitting, truncated probabilities, and label assignment."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy.special import logsumexp
from scipy.stats import truncnorm

from excavator2.fastcall import (
    assign_labels,
    correct_cellularity,
    fit_fastcall,
)

FIXTURES = Path(__file__).parent / "fixtures/legacy-fastcall"
FIT_CASES = [
    "uneven",
    "all_normal",
    "single",
    "nondefault",
    "paired_baseline",
]
PARAMETERS = {"rtol": 1e-13, "atol": 1e-14}
PROBABILITIES = {"rtol": 1e-13, "atol": 2e-15}


def load_case(name):
    path = (
        FIXTURES / "expected.npz" if name == "five_states" else FIXTURES / "edges" / f"{name}.npz"
    )
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def oracle_probabilities(values, fit):
    logs = np.full((len(values), 5), -np.inf)
    for j, (lower, upper) in enumerate(fit.bounds):
        if fit.priors[j] == 0:
            continue
        inside = (values >= lower) & (values <= upper)
        if np.isposinf(fit.deviations[j]):
            density = -np.log(upper - lower)
        else:
            sd, mean = fit.deviations[j], fit.means[j]
            density = truncnorm.logpdf(
                values[inside], (lower - mean) / sd, (upper - mean) / sd, loc=mean, scale=sd
            )
        logs[inside, j] = np.log(fit.priors[j]) + density
    totals = logsumexp(logs, axis=1)
    return np.exp(logs - totals[:, None]), totals.sum()


@pytest.mark.parametrize("name", ["five_states", *FIT_CASES])
@pytest.mark.parametrize("backend", ["python", "native"])
def test_fit_parameters_and_truncated_probabilities(name, backend):
    expected = load_case(name)
    fit = fit_fastcall(
        expected["mdata"],
        upper=expected["thru"][0],
        lower=expected["thrd"][0],
        backend=backend,
    )
    assert fit.converged
    assert_array_equal(fit.means, expected["muvec"])
    assert_array_equal(fit.bounds, expected["bound"])
    probabilities, likelihood = oracle_probabilities(expected["mdata"], fit)
    assert_allclose(fit.posterior, probabilities, **PROBABILITIES)
    assert_allclose(fit.trace[-1, -1], likelihood, rtol=1e-12, atol=1e-12)
    assert np.all(np.diff(fit.trace[:, -1]) >= -1e-10)
    assert_allclose(fit.priors.sum(), 1)


@pytest.mark.parametrize("name", ["boundaries", "extreme"])
@pytest.mark.parametrize("backend", ["python", "native"])
def test_out_of_support_fixtures_are_rejected(name, backend):
    expected = load_case(name)
    with pytest.raises(FloatingPointError, match="model support"):
        fit_fastcall(expected["mdata"], backend=backend)


@pytest.mark.parametrize("backend", ["python", "native"])
def test_corrected_fit_resolves_shared_nondefault_boundary(backend):
    case = load_case("nondefault")
    fit = fit_fastcall(case["mdata"], upper=case["thru"][0], lower=case["thrd"][0], backend=backend)
    index = np.flatnonzero(case["mdata"] == case["thru"][0])
    assert_array_equal(assign_labels(fit.posterior).labels[index], 1)
    assert np.all(fit.posterior[index, 3] > fit.posterior[index, 2])


def test_random_near_ties_replay_original_r_stream():
    expected = load_case("ties")
    result = assign_labels(expected["posterior"], r_seed=expected["seed_before"])
    assert_array_equal(result.labels, expected["calls"][:, 0])
    assert_array_equal(result.probabilities, expected["calls"][:, 1])
    assert_array_equal(result.r_seed, expected["seed_after"])
    # Streaming calls must consume the same words, including ties that later lose.
    first = assign_labels(expected["posterior"][:3], r_seed=expected["seed_before"])
    rest = assign_labels(expected["posterior"][3:], r_seed=first.r_seed)
    assert_array_equal(np.r_[first.labels, rest.labels], result.labels)
    assert_array_equal(rest.r_seed, result.r_seed)


def test_ties_require_explicit_historical_state():
    with pytest.raises(ValueError, match="original R RNG state"):
        assign_labels(load_case("ties")["posterior"])


def test_unique_calls_do_not_require_rng_state():
    expected = load_case("five_states")
    result = assign_labels(expected["posterior"])
    assert_array_equal(result.labels, expected["calls"][:, 0])
    assert result.r_seed is None


def test_cellularity_matches_original_and_preserves_input():
    expected = load_case("cellularity")
    values = expected["original"].copy()
    assert_allclose(
        correct_cellularity(values, expected["cellularity"][0]), expected["corrected"], **PARAMETERS
    )
    assert_array_equal(values, expected["original"])
    assert_array_equal(correct_cellularity(values), values)


@pytest.mark.parametrize("values", [[], [[0.0]], [np.nan], [np.inf]])
def test_invalid_segment_values_fail_explicitly(values):
    with pytest.raises(ValueError, match="segment values"):
        fit_fastcall(values)


def test_unsupported_rng_kind_fails_explicitly():
    state = load_case("ties")["seed_before"].copy()
    state[0] = 10407
    with pytest.raises(ValueError, match="Mersenne-Twister"):
        assign_labels(load_case("ties")["posterior"], r_seed=state)


@pytest.mark.parametrize("word", [-(2**31) - 1, 2**31])
def test_rng_state_rejects_out_of_range_words(word):
    state = load_case("ties")["seed_before"].astype(np.int64)
    state[2] = word
    with pytest.raises(ValueError, match="signed 32-bit"):
        assign_labels(load_case("ties")["posterior"], r_seed=state)


def test_golden_fixture_checksums():
    original = json.loads((FIXTURES / "manifest.json").read_text())
    assert (
        hashlib.sha256((FIXTURES / "expected.npz").read_bytes()).hexdigest()
        == original["expected_npz_sha256"]
    )
    edges = json.loads((FIXTURES / "edges/manifest.json").read_text())
    for name, record in edges["cases"].items():
        assert (
            hashlib.sha256((FIXTURES / "edges" / f"{name}.npz").read_bytes()).hexdigest()
            == record["sha256"]
        )


@pytest.mark.parametrize("backend", ["python", "native"])
@pytest.mark.parametrize("case", ["false_amplification", "missed_amplification"])
def test_a19_uses_fitted_support_for_final_calls(backend, case):
    if case == "false_amplification":
        values = np.r_[np.zeros(10), -0.4, 0.4, 1, 2, 3, 4, 5]
        labels = np.r_[np.zeros(11, dtype=int), 1, np.full(5, 2)]
    else:
        values = np.r_[np.tile([-0.3, 0.3], 50), 0.91, 5.0]
        labels = np.r_[np.zeros(100, dtype=int), 2, 2]
    fit = fit_fastcall(values, backend=backend)
    calls = assign_labels(fit.posterior)
    assert_array_equal(calls.labels, labels)
    assert_array_equal(calls.probabilities, 1)


@pytest.mark.parametrize("backend", ["python", "native"])
@pytest.mark.parametrize("upper, lower", [(0.35, 0.5), (0.2, 0.7)])
def test_truncated_posterior_matches_independent_density_at_final_parameters(backend, upper, lower):
    edges = np.array([-1.3, -lower, upper, 0.9])
    values = np.r_[
        -3.2,
        -3,
        -2.8,
        -1,
        -0.2,
        0,
        0.58,
        1,
        2,
        np.nextafter(edges, -np.inf),
        edges,
        np.nextafter(edges, np.inf),
    ]
    fit = fit_fastcall(values, backend=backend, upper=upper, lower=lower)
    expected, _ = oracle_probabilities(values, fit)
    inside = (values[:, None] >= fit.bounds[:, 0]) & (values[:, None] <= fit.bounds[:, 1])
    assert_allclose(fit.posterior, expected, rtol=1e-12, atol=2e-15)
    assert_array_equal(fit.posterior[~inside], 0)
    assert_allclose(fit.posterior.sum(axis=1), 1, rtol=0, atol=2e-15)


@pytest.mark.parametrize("backend", ["python", "native"])
def test_tiny_initial_width_stays_within_support(backend):
    fit = fit_fastcall([0.3, 0.3], backend=backend)
    assert_array_equal(assign_labels(fit.posterior).labels, [0, 0])
