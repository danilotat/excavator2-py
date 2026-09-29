"""Numerical contracts measured against the pinned, unmodified R implementation."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2.fastcall import (
    assign_labels,
    correct_cellularity,
    fit_fastcall,
    start_conditions,
    stopping_statistic,
)
from excavator2.reference.fastcall import posterior

FIXTURES = Path(__file__).parent / "fixtures/legacy-fastcall"
FIT_CASES = [
    "uneven",
    "boundaries",
    "all_normal",
    "single",
    "extreme",
    "nondefault",
    "paired_baseline",
]
# Observed max trace error is 7.1e-15, posterior error 2.3e-16 on macOS ARM64.
# Bounds below allow small library/reduction drift, never different iterations/calls.
PARAMETERS = {"rtol": 1e-13, "atol": 1e-14}
PROBABILITIES = {"rtol": 1e-13, "atol": 2e-15}


def load_case(name):
    path = (
        FIXTURES / "expected.npz" if name == "five_states" else FIXTURES / "edges" / f"{name}.npz"
    )
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


@pytest.mark.parametrize("name", ["five_states", *FIT_CASES])
@pytest.mark.parametrize("backend", ["python", "native"])
def test_fit_matches_original_r(name, backend):
    expected = load_case(name)
    fit = fit_fastcall(
        expected["mdata"], upper=expected["thru"][0], lower=expected["thrd"][0], backend=backend
    )
    assert fit.iterations == expected["iterations"][0]
    assert fit.converged
    assert_array_equal(fit.means, expected["muvec"])
    assert_array_equal(fit.bounds, expected["bound"])
    assert_allclose(fit.deviations, expected["sdvec"], **PARAMETERS)
    assert_allclose(fit.priors, expected["prior"], **PARAMETERS)
    assert_allclose(fit.posterior, expected["posterior"], **PROBABILITIES)
    if "iteration_trace" in expected:
        assert_allclose(fit.trace, expected["iteration_trace"], **PARAMETERS)
        means, deviations = start_conditions(
            expected["mdata"], expected["thru"][0], expected["thrd"][0]
        )
        assert_allclose(deviations, expected["initial_deviations"], **PARAMETERS)
        priors = np.array([0.05, 0.1, 0.7, 0.1, 0.05])
        statistic = stopping_statistic(
            posterior(expected["mdata"], means, deviations, priors), priors
        )
        assert_allclose(statistic, expected["initial_statistic"][0], **PARAMETERS)
    assigned = assign_labels(fit.posterior, r_seed=expected["seed_before"])
    assert_array_equal(assigned.labels, expected["calls"][:, 0])
    assert_allclose(assigned.probabilities, expected["calls"][:, 1], **PROBABILITIES)
    assert_array_equal(assigned.r_seed, expected["seed_after"])


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
