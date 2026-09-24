"""Python FastCall policy and iteration loop, preserving the original R behavior.

This API accepts already constructed segment values. Segmentation, segment-table
construction, CN/VCF rendering, and CLI integration are separate future stages.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .reference import fastcall as kernels

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class FastCallFit:
    means: FloatArray
    deviations: FloatArray
    priors: FloatArray
    bounds: FloatArray
    iterations: int
    converged: bool
    posterior: FloatArray
    # Each row: five deviations, five priors, legacy stopping statistic.
    trace: FloatArray


@dataclass(frozen=True)
class FastCallLabels:
    labels: NDArray[np.int32]
    probabilities: FloatArray
    # Updated R MT state, only when a replay state was provided.
    r_seed: NDArray[np.int32] | None


def segment_values(values: ArrayLike) -> FloatArray:
    """Validate the finite, nonempty segment domain supported by this reference."""
    result = np.asarray(values, dtype=np.float64)
    if result.ndim != 1 or result.size == 0 or not np.isfinite(result).all():
        raise ValueError("segment values must be a nonempty, finite one-dimensional array")
    return result


def correct_cellularity(values: ArrayLike, cellularity: float = 1.0) -> FloatArray:
    """Apply the original correction/floor only to the calling values."""
    values = segment_values(values)
    if not np.isfinite(cellularity) or not 0 < cellularity <= 1:
        raise ValueError("cellularity must be in (0, 1]")
    if cellularity == 1:
        return values.copy()
    with np.errstate(over="ignore"):
        corrected = np.exp2(values) / cellularity - (1 - cellularity) / cellularity
    corrected[corrected < 2**-5] = 2**-5
    return segment_values(np.log2(corrected))


def start_conditions(
    values: FloatArray, upper: float, lower: float
) -> tuple[FloatArray, FloatArray]:
    """StartCond: fixed means, sample SD within overlapping inclusive ranges."""
    means = np.array([-3.0, -1.0, 0.0, 0.58, 1.0])
    deviations = np.full(5, 0.01)
    boundaries = [-50.0, -1.5, -lower, upper, 0.9, 50.0]
    for j in range(5):
        selected = values[(values >= boundaries[j]) & (values <= boundaries[j + 1])]
        if selected.size > 1:
            deviations[j] = np.std(selected, ddof=1)
    deviations[deviations < 0.001] = 0.001
    return means, deviations


def stopping_statistic(posterior: FloatArray, priors: FloatArray) -> float:
    """R sum(PosteriorP(...) * prior), including column-major vector recycling.

    This is intentionally NOT a log likelihood or a column-weighted sum.
    """
    flat = posterior.ravel(order="F")
    recycled = np.resize(priors, flat.size)
    return float(np.sum(flat * recycled))


def fit_fastcall(values: ArrayLike, *, upper: float = 0.35, lower: float = 0.5) -> FastCallFit:
    """Fit the legacy five-state model to one value per segment.

    `lower` is the positive magnitude d in the legacy YAML; the normal state's
    lower bound is -d. No segment-length weighting or learned means are added.
    """
    values = segment_values(values)
    if not np.isfinite([upper, lower]).all() or not (0 < upper < 0.9 and 0 < lower < 1.3):
        raise ValueError("require 0 < upper < 0.9 and 0 < lower < 1.3")
    means, deviations = start_conditions(values, upper, lower)
    priors = np.array([0.05, 0.1, 0.7, 0.1, 0.05])
    edges = [-20.0, -1.3, -lower, upper, 0.9, 20.0]
    bounds = np.column_stack([edges[:-1], edges[1:]])
    probabilities = kernels.posterior(values, means, deviations, priors)
    statistic = stopping_statistic(probabilities, priors)
    trace = []
    converged = False
    for iteration in range(1, 1001):
        responsibilities = kernels.expectation(values, means, deviations, priors, bounds)
        deviations, priors = kernels.maximization(values, responsibilities, means, deviations)
        probabilities = kernels.posterior(values, means, deviations, priors)
        previous = statistic
        statistic = stopping_statistic(probabilities, priors)
        if not np.isfinite(statistic):
            raise FloatingPointError("legacy FastCall produced a non-finite stopping statistic")
        trace.append(np.r_[deviations, priors, statistic])
        if abs(statistic - previous) < 1e-5:
            converged = True
            break
    return FastCallFit(
        means, deviations, priors, bounds, iteration, converged, probabilities, np.array(trace)
    )


class _RUniformReplay:
    """Use NumPy's MT engine with an explicit R .Random.seed, not NumPy seeding.

    R consumes one 32-bit MT word per uniform. NumPy's usual floating-point
    generator consumes a different number of words, so it cannot be used here.
    """

    def __init__(self, seed: ArrayLike):
        raw = np.asarray(seed)
        if raw.shape != (626,) or raw.dtype.kind not in "iu":
            raise ValueError("expected the 626 integers of an R Mersenne-Twister .Random.seed")
        self.kind = int(raw[0])
        position = int(raw[1])
        if self.kind % 100 != 3 or not 0 <= position <= 624:
            raise ValueError("only an R Mersenne-Twister state is supported")
        self.engine = np.random.RandomState()
        self.engine.set_state(("MT19937", raw[2:].astype(np.uint32), position, 0, 0.0))

    def uniform(self) -> float:
        word = int(self.engine.randint(0, 2**32, dtype=np.uint32))
        # R excludes the uniform endpoints; 1 is impossible for a uint32 word.
        return word / 2**32 if word else 0.5 / (2**32 - 1)

    def state(self) -> NDArray[np.int32]:
        _, words, position, _, _ = self.engine.get_state()
        return np.r_[np.array([self.kind, position], dtype=np.int32), words.view(np.int32)]


def assign_labels(posterior: ArrayLike, *, r_seed: ArrayLike | None = None) -> FastCallLabels:
    """Assign legacy labels, replaying R's random near-tie policy when needed.

    With no R state, only unambiguous maxima are accepted. A near tie raises
    instead of silently choosing a different call or inventing a historical seed.
    With a state, even intermediate losing ties consume the original uniforms.
    """
    probabilities = np.asarray(posterior, dtype=np.float64)
    if (
        probabilities.ndim != 2
        or probabilities.shape[1] != 5
        or not np.isfinite(probabilities).all()
        or (probabilities < 0).any()
    ):
        raise ValueError("posterior must be a finite nonnegative (n, 5) matrix")
    replay = _RUniformReplay(r_seed) if r_seed is not None else None
    indices = np.empty(probabilities.shape[0], dtype=np.int32)
    for i, row in enumerate(probabilities):
        tolerance = 1e-5 * np.max(row)
        if replay is None:
            top = np.sort(row)[-2:]
            if top[1] - top[0] <= tolerance:
                raise ValueError(
                    f"row {i} has a legacy random tie; provide its original R RNG state"
                )
            indices[i] = np.argmax(row)
            continue
        best_value, best_index, tied = row[0], 0, 1
        for j in range(1, 5):
            if row[j] > best_value + tolerance:
                best_value, best_index, tied = row[j], j, 1
            elif row[j] >= best_value - tolerance:
                tied += 1
                if tied * replay.uniform() < 1:
                    best_index = j
        indices[i] = best_index
    return FastCallLabels(
        indices - 2,
        probabilities[np.arange(len(indices)), indices],
        replay.state() if replay is not None else None,
    )
