"""Differential checks against the pinned original R/Fortran implementation."""

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2 import _core
from excavator2.hslm import (
    distance_covariates,
    estimate_parameters,
    filter_breaks,
    reconstruct,
    segment,
    state_grid,
)
from excavator2.reference import hslm as reference

FIXTURES = Path(__file__).parent / "fixtures" / "legacy-hslm"


@pytest.mark.parametrize("case", ["shifts", "short", "uneven", "single", "ties"])
def test_original_kernel_traces(case):
    with np.load(FIXTURES / f"{case}.npz") as f:

        def scalar(name):
            return float(f[name][0])

        values, means, eta, initial = (
            np.ascontiguousarray(f[name]) for name in ["values", "means", "eta", "initial"]
        )
        args = (values, means, scalar("mi"), scalar("smu"), scalar("sepsilon"), eta)
        transitions, emissions = reference.matrices(*args)
        path, scores, predecessors = reference.viterbi(initial, transitions, emissions)
        k, n = len(means), len(values)
        expected_transitions = f["transitions"].T.reshape(n - 1, k, k).transpose(0, 2, 1)
        for actual, expected in [
            (transitions, expected_transitions),
            (emissions, f["emissions"].T),
            (scores, f["scores"].T),
        ]:
            assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)
        assert_array_equal(path, f["path"])
        assert_array_equal(predecessors, f["predecessors"].T)
        for trace in [False, True]:
            native = _core.hslm_segment(*args, initial, trace)
            assert_array_equal(native[0], f["path"])
            if trace:
                for actual, expected in zip(
                    native[1:4], [transitions, emissions, scores], strict=True
                ):
                    assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)
                assert_array_equal(native[4], predecessors)
        assert_array_equal(filter_breaks(f["breaks"], scalar("fw")), f["filtered"])
        assert_array_equal(reconstruct(values, f["filtered"].astype(int)), f["segmented"])
        estimated = estimate_parameters(f["estimation_values"], scalar("omega"))
        assert_allclose(
            [estimated.mi, estimated.smu, estimated.sepsilon],
            [scalar("mi"), scalar("smu"), scalar("sepsilon")],
            rtol=2e-14,
        )
        assert_allclose(
            distance_covariates(f["positions"], scalar("theta"), scalar("distance")),
            eta,
            rtol=2e-14,
        )
        if case not in ("ties", "single"):
            assert_array_equal(state_grid(), means)
            for backend in ["python", "native"]:
                result = segment(
                    values,
                    f["positions"],
                    estimated,
                    theta=scalar("theta"),
                    distance=scalar("distance"),
                    min_windows=scalar("fw"),
                    backend=backend,
                    trace=True,
                )
                assert_array_equal(result.path, path)
                assert_array_equal(result.values, f["segmented"])


def test_invalid_domains_fail_explicitly():
    with pytest.raises(ValueError, match="variance"):
        estimate_parameters(np.zeros(10))
    with pytest.raises(ValueError, match="nondecreasing"):
        distance_covariates([2, 1], 0.1, 100)
    p = estimate_parameters(np.linspace(-1, 1, 10))
    with pytest.raises(ValueError, match="single-window"):
        segment([0], [1], p)
    with pytest.raises(ValueError, match="same length"):
        segment([0, 1], [1], p)


@pytest.mark.parametrize(
    "case", ["dtype", "stride", "alignment", "empty", "rank", "eta", "sd", "initial", "nan"]
)
def test_native_rejects_invalid_buffers(case):
    values = np.zeros(4)
    means = np.array([-0.1, 0.1])
    eta = np.full(3, 0.1)
    initial = np.full(2, np.log(0.5))
    sd = 0.1
    exception = ValueError
    if case == "dtype":
        values = values.astype(np.float32)
        exception = TypeError
    elif case == "stride":
        values = np.zeros(8)[::2]
        exception = TypeError
    elif case == "alignment":
        values = np.ndarray((4,), dtype=np.float64, buffer=bytearray(33), offset=1)
    elif case == "empty":
        values = np.array([])
    elif case == "rank":
        values = values.reshape(2, 2)
    elif case == "eta":
        eta[1] = 0
    elif case == "sd":
        sd = -1
    elif case == "initial":
        initial = initial[:1]
    else:
        values[0] = np.nan
    with pytest.raises(exception):
        _core.hslm_segment(values, means, 0, sd, sd, eta, initial)
