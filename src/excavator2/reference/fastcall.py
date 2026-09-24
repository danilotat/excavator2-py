"""Readable numerical kernels matching LibraryFastCall.R.

Inputs are binary64 arrays: segment values (n,), parameters (5,), and
responsibilities (n, 5). Bounds are inclusive. These functions deliberately
retain the legacy zero-density fallback and fixed component means.
"""

import numpy as np
from numpy.typing import NDArray
from scipy.special import ndtr

FloatArray = NDArray[np.float64]


def normal_density(values: FloatArray, mean: float, deviation: float) -> FloatArray:
    """Normal density in ordinary probability space, as used by the original."""
    standardized = (values - mean) / deviation
    return np.exp(-0.5 * standardized**2) / np.sqrt(2.0 * np.pi) / deviation


def normalize_weights(weights: FloatArray, values: FloatArray, means: FloatArray) -> FloatArray:
    """Normalize rows; underflow chooses the first closest fixed mean.

    R rowSums(na.rm=TRUE) ignores missing terms but leaves them in the numerator.
    Do not substitute a log-space normalization: that changes legacy underflow.
    """
    totals = np.nansum(weights, axis=1)
    empty = np.flatnonzero(totals == 0)
    if empty.size:
        nearest = np.argmin(np.abs(values[empty, None] - means), axis=1)
        weights[empty, nearest] = 1.0
        totals = np.nansum(weights, axis=1)
    return weights / totals[:, None]


def posterior(
    values: FloatArray, means: FloatArray, deviations: FloatArray, priors: FloatArray
) -> FloatArray:
    """Original PosteriorP: final densities are not truncated."""
    weights = np.column_stack(
        [priors[j] * normal_density(values, means[j], deviations[j]) for j in range(5)]
    )
    return normalize_weights(weights, values, means)


def expectation(
    values: FloatArray,
    means: FloatArray,
    deviations: FloatArray,
    priors: FloatArray,
    bounds: FloatArray,
) -> FloatArray:
    """Original EStep/gfct: truncated densities, inclusive bounds, Inf -> 100."""
    densities = np.zeros((values.size, 5), dtype=np.float64)
    for j, (lower, upper) in enumerate(bounds):
        inside = (values >= lower) & (values <= upper)
        if not inside.any():
            continue
        denominator = ndtr((upper - means[j]) / deviations[j]) - ndtr(
            (lower - means[j]) / deviations[j]
        )
        with np.errstate(divide="ignore", invalid="ignore", under="ignore"):
            column = normal_density(values, means[j], deviations[j]) * inside / denominator
        column[np.isposinf(column)] = 100.0
        densities[:, j] = column
    return normalize_weights(densities * priors, values, means)


def maximization(
    values: FloatArray, responsibilities: FloatArray, means: FloatArray, deviations: FloatArray
) -> tuple[FloatArray, FloatArray]:
    """Original MStep updates deviations/priors only; means stay fixed."""
    mass = np.sum(responsibilities, axis=0)
    new_deviations = deviations.copy()
    new_priors = np.empty(5, dtype=np.float64)
    for j in range(5):
        if mass[j] != 0:
            deviation = np.sqrt(np.dot(responsibilities[:, j], (values - means[j]) ** 2) / mass[j])
            if not deviation < 1e-100:
                new_deviations[j] = deviation
            new_priors[j] = mass[j] / values.size
        else:
            new_priors[j] = 1e-6
    return new_deviations, new_priors / np.sum(new_priors)
