"""Fit and classify segment values with truncated FastCall components."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import ndtr

from . import _core
from .reference import fastcall as reference_kernels

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class FastCallFit:
    means: FloatArray
    # Positive infinity denotes a uniform density within the class bounds.
    deviations: FloatArray
    priors: FloatArray
    bounds: FloatArray
    iterations: int
    converged: bool
    posterior: FloatArray
    # Each row: five deviations, five priors, truncated log-likelihood.
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
    """Correct segment ratios for both calling and diploid-equivalent reporting."""
    values = segment_values(values)
    if not np.isfinite(cellularity) or not 0 < cellularity <= 1:
        raise ValueError("cellularity must be in (0, 1]")
    if cellularity == 1:
        return values.copy()
    with np.errstate(over="ignore"):
        corrected = np.exp2(values) / cellularity - (1 - cellularity) / cellularity
    corrected[corrected < 2**-5] = 2**-5
    return segment_values(np.log2(corrected))


def recenter_ploidy(values: ArrayLike, ploidy: float = 2.0, cellularity: float = 1.0) -> FloatArray:
    """Shift median-centered mixture ratios so zero is diploid, given tumor ploidy."""
    values = segment_values(values)
    if not np.isfinite(ploidy) or ploidy <= 0:
        raise ValueError("ploidy must be positive and finite")
    if not np.isfinite(cellularity) or not 0 < cellularity <= 1:
        raise ValueError("cellularity must be in (0, 1]")
    if ploidy == 2:
        return values.copy()
    mixture = cellularity * ploidy + 2 * (1 - cellularity)
    return segment_values(values + np.log2(mixture / 2))


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


def class_edges(upper: float = 0.35, lower: float = 0.5) -> FloatArray:
    """Interior class boundaries on the cellularity-corrected log2 scale."""
    if not np.isfinite([upper, lower]).all() or not (0 < upper < 0.9 and 0 < lower < 1.3):
        raise ValueError("require 0 < upper < 0.9 and 0 < lower < 1.3")
    return np.array([-1.3, -lower, upper, 0.9])


def window_noise(values: ArrayLike, classes: ArrayLike, starts, ends) -> dict[str, float]:
    """Robust per-window noise for IN and OUT from adjacent same-class differences.

    Only pairs inside one segment are used, so copy-number steps do not inflate
    the estimate. A class without such pairs borrows the estimate over all pairs.
    """
    values, classes = segment_values(values), np.asarray(classes)
    if classes.shape != values.shape or not np.isin(classes, ["IN", "OUT"]).all():
        raise ValueError("classes must contain IN/OUT labels matching values")
    segment = np.repeat(np.arange(len(starts)), np.subtract(ends, starts))
    if segment.shape != values.shape:
        raise ValueError("segments must partition the windows")

    def estimate(selected):
        index = np.flatnonzero(selected)
        same = segment[index[1:]] == segment[index[:-1]]
        differences = np.abs(np.diff(values[index]))[same]
        if len(differences) < 2:
            return None
        # Median |difference| of two independent N(0, sigma^2) windows is 0.9539 sigma.
        return float(max(np.median(differences) / 0.9538725524089398, 0.001))

    pooled = estimate(np.ones(len(values), dtype=bool))
    if pooled is None:
        raise ValueError("measurement noise requires adjacent windows within segments")
    return {kind: estimate(classes == kind) or pooled for kind in ["IN", "OUT"]}


def segment_standard_errors(classes: ArrayLike, starts, ends, noise) -> FloatArray:
    """Large-sample standard error of each segment median of independent windows."""
    classes = np.asarray(classes)
    starts, ends = np.asarray(starts), np.asarray(ends)
    variance = np.where(classes == "IN", noise["IN"] ** 2, noise["OUT"] ** 2)
    total = np.add.reduceat(variance, starts)
    return np.sqrt(np.pi / 2) * np.sqrt(total) / (ends - starts)


def call_probabilities(
    segments: ArrayLike,
    standard_errors: ArrayLike,
    *,
    upper: float = 0.35,
    lower: float = 0.5,
    cellularity: float = 1.0,
) -> FloatArray:
    """Probability that each true segment level lies in each class interval.

    The observed mixture median is normal around the true level with its standard
    error, under a flat prior. Corrected class edges map back to the mixture scale
    through the monotone cellularity correction; the outer classes are open-ended,
    so each row sums to one and mixtures below the correction floor are deep losses.
    """
    segments = segment_values(segments)
    errors = np.asarray(standard_errors, dtype=np.float64)
    if errors.shape != segments.shape or not np.isfinite(errors).all() or (errors <= 0).any():
        raise ValueError("standard errors must be positive, finite and match segments")
    if not np.isfinite(cellularity) or not 0 < cellularity <= 1:
        raise ValueError("cellularity must be in (0, 1]")
    edges = np.r_[-np.inf, class_edges(upper, lower), np.inf]
    if cellularity < 1:
        with np.errstate(divide="ignore"):
            edges = np.log2(cellularity * np.exp2(edges) + (1 - cellularity))
        edges[0] = -np.inf
    cumulative = ndtr((edges - segments[:, None]) / errors[:, None])
    return np.diff(cumulative, axis=1)


def fit_fastcall(
    values: ArrayLike,
    *,
    upper: float = 0.35,
    lower: float = 0.5,
    backend: str = "native",
) -> FastCallFit:
    """Fit five states; the normal interval is [-lower, upper]."""
    values = np.require(segment_values(values), dtype=np.float64, requirements=["C", "A"])
    if backend == "native":
        expectation = _core.fastcall_expectation
        maximization = _core.fastcall_maximization
    elif backend == "python":
        expectation = reference_kernels.expectation
        maximization = reference_kernels.maximization
    else:
        raise ValueError("backend must be 'native' or 'python'")
    if not np.isfinite([upper, lower]).all() or not (0 < upper < 0.9 and 0 < lower < 1.3):
        raise ValueError("require 0 < upper < 0.9 and 0 < lower < 1.3")
    means, deviations = start_conditions(values, upper, lower)
    priors = np.array([0.05, 0.1, 0.7, 0.1, 0.05])
    edges = np.r_[-20.0, class_edges(upper, lower), 20.0]
    bounds = np.column_stack([edges[:-1], edges[1:]])
    probabilities, likelihood = expectation(values, means, deviations, priors, bounds)
    trace = []
    converged = False
    for iteration in range(1, 1001):
        old_precision, old_priors, previous = 1 / deviations**2, priors, likelihood
        deviations, priors = maximization(values, probabilities, means, deviations, bounds)
        probabilities, likelihood = expectation(values, means, deviations, priors, bounds)
        if not np.isfinite(likelihood):
            raise FloatingPointError("non-finite truncated log-likelihood")
        if likelihood < previous - 1e-10 * (1 + abs(previous)):
            raise FloatingPointError("truncated log-likelihood decreased")
        trace.append(np.r_[deviations, priors, likelihood])
        precision_change = np.max(np.abs(1 / deviations**2 - old_precision) / (1 + old_precision))
        if (
            abs(likelihood - previous) <= 1e-8 * (1 + abs(previous))
            and precision_change <= 1e-8
            and np.max(np.abs(priors - old_priors)) <= 1e-8
        ):
            converged = True
            break
    return FastCallFit(
        means,
        deviations,
        priors,
        bounds,
        iteration,
        converged,
        probabilities,
        np.array(trace),
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
        if (raw < np.iinfo(np.int32).min).any() or (raw > np.iinfo(np.int32).max).any():
            raise ValueError("R RNG state entries must be signed 32-bit integers")
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
