"""Single-profile HSLM setup, segmentation and minimum-support filtering."""

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
    deviations: NDArray[np.float64] | None = None


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


def estimate_parameters(values, omega=0.1, classes=None):
    """Estimate arm noise from robust adjacent differences."""
    values = vector(values, "values")
    if not 0 < omega < 1:
        raise ValueError("omega must lie strictly between zero and one")

    def noise_estimate(observations):
        differences = np.abs(np.diff(observations))
        noise = np.median(differences) / 0.9538725524089398 if len(differences) else 0.0
        return float(max(noise, 0.001))

    noise = noise_estimate(values)
    deviations = None
    if classes is not None:
        classes = np.asarray(classes)
        if classes.shape != values.shape or not np.isin(classes, ["IN", "OUT"]).all():
            raise ValueError("classes must contain IN/OUT labels matching values")
        deviations = np.full(len(values), noise)
        for label in np.unique(classes):
            selected = classes == label
            if selected.sum() >= 3:
                deviations[selected] = noise_estimate(values[selected])
        noise = float(deviations.min())
    return Parameters(0.0, float(np.sqrt(omega)), noise, deviations)


def state_grid(values, noise):
    """Use occupied noise-scaled bins, with zero and no amplitude clipping."""
    if not np.isfinite(noise) or noise <= 0:
        raise ValueError("noise must be positive and finite")
    step = max(0.001, min(0.05, noise / 2))
    return np.unique(np.r_[0.0, np.round(vector(values, "values") / step) * step])


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


def filter_breaks(breaks, min_windows, values, support=None):
    """Merge insufficient-support segments toward the closest median."""
    values = vector(values, "values")
    if isinstance(min_windows, bool) or int(min_windows) != min_windows or min_windows < 0:
        raise ValueError("min_windows must be a nonnegative integer")
    result = np.asarray(breaks, dtype=np.int64).copy()
    if (
        result.ndim != 1
        or len(result) < 2
        or result[0] != 0
        or result[-1] != len(values)
        or (np.diff(result) <= 0).any()
    ):
        raise ValueError("breaks must partition all values")
    support = np.ones(len(values), dtype=bool) if support is None else np.asarray(support)
    if support.shape != values.shape or support.dtype != np.bool_:
        raise ValueError("support must be a boolean vector matching values")
    while len(result) > 2:
        counts = np.array([support[a:b].sum() for a, b in zip(result[:-1], result[1:])])
        short = counts < min_windows
        eligible = short[:-1] | short[1:]
        if not eligible.any():
            break
        medians = np.array([np.median(values[a:b]) for a, b in zip(result[:-1], result[1:])])
        costs = np.abs(np.diff(medians))
        # Remove all equally good boundaries together, including symmetric ties.
        chosen = eligible & (costs == costs[eligible].min())
        result = np.delete(result, np.flatnonzero(chosen) + 1)
    return result


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
    support=None,
    backend="native",
    trace=False,
):
    """Segment one arm with an explicit stationary boundary prior."""
    values = vector(values, "values")
    positions = vector(positions, "positions")
    if len(positions) != len(values):
        raise ValueError("values and positions must have the same length")
    if isinstance(min_windows, bool) or int(min_windows) != min_windows or min_windows < 0:
        raise ValueError("min_windows must be a nonnegative integer")
    p = parameters
    if not np.isfinite([p.mi, p.smu, p.sepsilon]).all() or p.smu <= 0 or p.sepsilon <= 0:
        raise ValueError("parameters require finite means and positive deviations")
    deviations = None
    if p.deviations is not None:
        deviations = vector(p.deviations, "deviations")
        if len(deviations) != len(values) or (deviations <= 0).any():
            raise ValueError("deviations must be positive and match values")
    eta = distance_covariates(positions, theta, distance)
    means = state_grid(values, p.sepsilon)
    initial = -((means - p.mi) ** 2 / (2 * p.smu**2))
    norm = float(initial[0])
    for value in initial[1:]:
        norm = reference.elnsum(norm, float(value))
    initial -= norm
    if backend == "python":
        transitions, emissions = reference.matrices(
            values, means, p.mi, p.smu, p.sepsilon, eta, deviations
        )
        path, scores, predecessors = reference.viterbi(initial, transitions, emissions)
    elif backend == "native":
        path, transitions, emissions, scores, predecessors = _core.hslm_segment(
            values, means, p.mi, p.smu, p.sepsilon, eta, initial, trace, deviations
        )
    else:
        raise ValueError("backend must be 'native' or 'python'")
    breaks = np.r_[0, np.flatnonzero(np.diff(path)) + 1, len(path)]
    filtered = filter_breaks(breaks, min_windows, values, support)
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
