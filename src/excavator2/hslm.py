"""Readable single-profile HSLM setup, segmentation and legacy filtering policy."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from . import _core
from .reference import hslm as reference


@dataclass(frozen=True)
class Parameters:
    mi: float
    smu: float
    sepsilon: float


@dataclass(frozen=True)
class Segmentation:
    path: NDArray[np.int32]
    breaks: NDArray[np.int64]
    filtered_breaks: NDArray[np.int64]
    values: NDArray[np.float64]
    transitions: NDArray[np.float64] | None = None
    emissions: NDArray[np.float64] | None = None
    scores: NDArray[np.float64] | None = None
    predecessors: NDArray[np.int32] | None = None


def vector(values, name):
    result = np.require(values, dtype=np.float64, requirements=["C", "A"])
    if result.ndim != 1 or not len(result) or not np.isfinite(result).all():
        raise ValueError(f"{name} must be a finite nonempty vector")
    return result


def estimate_parameters(values, omega=0.1):
    """R ParamEstSeq: global 1–99% inclusive trimming and sample variance."""
    values = vector(values, "values")
    if not 0 < omega < 1:
        raise ValueError("omega must lie strictly between zero and one")
    low, high = np.quantile(values, [0.01, 0.99], method="linear")
    selected = values[(values >= low) & (values <= high)]
    if len(selected) < 2:
        raise ValueError("legacy parameter estimation requires two retained values")
    variance = np.var(selected, ddof=1)
    if not np.isfinite(variance) or variance <= 0:
        raise ValueError("legacy HSLM requires positive finite global variance")
    return Parameters(0.0, float(np.sqrt(omega * variance)), float(np.sqrt((1 - omega) * variance)))


def state_grid():
    """R seq(-1, 1, by=0.1), including the binary64 operation order."""
    return -1.0 + np.arange(21, dtype=np.float64) * 0.1


def distance_covariates(positions, theta, distance):
    positions = vector(positions, "positions")
    if not 0 < theta < 1 or not np.isfinite(distance) or distance <= 0:
        raise ValueError("require 0 < theta < 1 and a positive finite distance")
    # The inference wrapper uses as.integer, truncating toward zero before diff.
    positions = np.trunc(positions)
    if (np.abs(positions) > np.iinfo(np.int32).max).any():
        raise ValueError("positions exceed the legacy integer range")
    delta = np.diff(positions)
    if (delta < 0).any():
        raise ValueError("positions must be nondecreasing within an arm")
    with np.errstate(divide="ignore"):
        return theta + (1 - theta) * np.exp(np.log(theta) / (delta / distance))


def filter_breaks(breaks, min_windows):
    """FilterSeg removes starts of short segments, with its first-segment quirk.

    If the sole segment is short, the final endpoint is removed and reconstruction
    remains zero. Preserve this surprising legacy behavior until compatibility.
    """
    result = np.asarray(breaks, dtype=np.int64)
    short = np.flatnonzero(np.diff(result) <= min_windows)
    if len(short) and short[0] == 0:
        short[0] = 1
        short = np.unique(short)
    return np.delete(result, short)


def reconstruct(values, breaks):
    result = np.zeros(len(values))
    for start, end in zip(breaks[:-1], breaks[1:], strict=True):
        result[start:end] = np.median(values[start:end])
    return result


def segment(
    values,
    positions,
    parameters,
    *,
    theta=1e-5,
    distance=1e6,
    min_windows=2,
    backend="native",
    trace=False,
):
    """Segment one arm using parameters estimated over the complete profile.

    Trace matrices are optional; the native recurrence uses rolling scores and
    on-demand transition blocks for ordinary calls. A single-window arm is rejected
    because legacy SortState fails on it; no new scientific behavior is substituted.
    """
    values = vector(values, "values")
    if len(values) < 2:
        raise ValueError("legacy SortState does not support a single-window arm")
    positions = vector(positions, "positions")
    if len(positions) != len(values):
        raise ValueError("values and positions must have the same length")
    if isinstance(min_windows, bool) or int(min_windows) != min_windows or min_windows < 0:
        raise ValueError("min_windows must be a nonnegative integer")
    p = parameters
    if not np.isfinite([p.mi, p.smu, p.sepsilon]).all() or p.smu <= 0 or p.sepsilon <= 0:
        raise ValueError("parameters require finite means and positive deviations")
    eta = distance_covariates(positions, theta, distance)
    means = state_grid()
    initial = np.full(len(means), np.log(1 / len(means)))
    if backend == "python":
        transitions, emissions = reference.matrices(values, means, p.mi, p.smu, p.sepsilon, eta)
        path, scores, predecessors = reference.viterbi(initial, transitions, emissions)
    elif backend == "native":
        path, transitions, emissions, scores, predecessors = _core.hslm_segment(
            values, means, p.mi, p.smu, p.sepsilon, eta, initial, trace
        )
    else:
        raise ValueError("backend must be 'native' or 'python'")
    breaks = np.r_[0, np.flatnonzero(np.diff(path)) + 1, len(path)]
    filtered = filter_breaks(breaks, min_windows)
    return Segmentation(
        path,
        breaks,
        filtered,
        reconstruct(values, filtered),
        transitions if trace else None,
        emissions if trace else None,
        scores if trace else None,
        predecessors if trace else None,
    )
