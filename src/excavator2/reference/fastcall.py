"""Log-space probabilities and fixed-mean truncated-normal fitting."""

import numpy as np
from numpy.typing import NDArray
from scipy.special import erf, erfcx, logsumexp

FloatArray = NDArray[np.float64]
MIN_DEVIATION = 0.001
LOG_ROOT_2PI = 0.5 * np.log(2 * np.pi)
NODES, WEIGHTS = np.polynomial.legendre.leggauss(16)


def log_scaled_mass(a: float, b: float) -> float:
    """Log integral relative to the largest density within the interval."""
    if not np.isfinite([a, b]).all() or not a < b:
        raise FloatingPointError("unrepresentable truncated interval")
    if b <= 0:
        a, b = -b, -a
    width = b - a
    if width * max(1.0, abs(a), abs(b)) < 1:
        anchor = max(0.0, a)
        offsets = width / 2 * (NODES + 1) + (a - anchor)
        return float(
            np.log(width / 2)
            + np.log(np.dot(WEIGHTS, np.exp(-anchor * offsets - 0.5 * offsets**2)))
        )
    if a < 0:
        return float(np.log(np.sqrt(np.pi / 2) * (erf(b / np.sqrt(2)) - erf(a / np.sqrt(2)))))
    left = np.log(np.sqrt(np.pi / 2) * erfcx(a / np.sqrt(2)))
    right = np.log(np.sqrt(np.pi / 2) * erfcx(b / np.sqrt(2)))
    ratio = -0.5 * (b - a) * (b + a) + right - left
    return float(left + np.log(-np.expm1(ratio)))


def normalize_logs(log_weights: FloatArray) -> tuple[FloatArray, float]:
    totals = logsumexp(log_weights, axis=1)
    if not np.isfinite(totals).all():
        raise FloatingPointError("observation has no finite probability within model support")
    # Subtract the maximum first to retain accuracy for very negative logs.
    shifted = log_weights - np.max(log_weights, axis=1, keepdims=True)
    weights = np.exp(shifted)
    with np.errstate(over="ignore"):
        likelihood = float(np.sum(totals))
    if not np.isfinite(likelihood):
        raise FloatingPointError("non-finite mixture log-likelihood")
    return weights / weights.sum(axis=1, keepdims=True), likelihood


def posterior(
    values: FloatArray, means: FloatArray, deviations: FloatArray, priors: FloatArray
) -> FloatArray:
    """Untruncated Gaussian probabilities, normalized in log space."""
    if not np.isfinite(deviations).all():
        raise ValueError("untruncated deviations must be finite")
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        weights = (
            np.log(priors)
            - np.log(deviations)
            - LOG_ROOT_2PI
            - 0.5 * ((values[:, None] - means) / deviations) ** 2
        )
    return normalize_logs(weights)[0]


def expectation(
    values: FloatArray,
    means: FloatArray,
    deviations: FloatArray,
    priors: FloatArray,
    bounds: FloatArray,
) -> tuple[FloatArray, float]:
    """Return responsibilities and truncated mixture log-likelihood."""
    logs = np.full((values.size, 5), -np.inf)
    for j, (lower, upper) in enumerate(bounds):
        if not lower < upper:
            raise ValueError("bounds must have positive width")
        inside = (values >= lower) & (values <= upper)
        if priors[j] == 0 or not inside.any():
            continue
        sd = deviations[j]
        if np.isposinf(sd):
            logs[inside, j] = np.log(priors[j]) - np.log(upper - lower)
        else:
            mass = log_scaled_mass((lower - means[j]) / sd, (upper - means[j]) / sd)
            anchor = np.clip(means[j], lower, upper)
            difference = (values[inside] - anchor) / sd
            with np.errstate(over="ignore", invalid="ignore"):
                logs[inside, j] = (
                    np.log(priors[j])
                    - np.log(sd)
                    - mass
                    - 0.5 * difference * (difference + 2 * (anchor - means[j]) / sd)
                )
    return normalize_logs(logs)


def second_moment(sd: float, lower: float, upper: float) -> float:
    """E[(X-mu)^2] for bounds already centered on mu."""
    if np.isposinf(sd):
        return (lower**2 + lower * upper + upper**2) / 3
    a, b = lower / sd, upper / sd
    if max(abs(a), abs(b)) <= 1 or (b - a) * max(1.0, abs(a), abs(b)) < 1:
        mid, half = lower + (upper - lower) / 2, (upper - lower) / 2
        offsets = half * NODES
        weights = WEIGHTS * np.exp(-(mid / sd) * (offsets / sd) - 0.5 * (offsets / sd) ** 2)
        return float(np.dot(weights, (mid + offsets) ** 2) / weights.sum())
    mass = log_scaled_mass(a, b)
    anchor = np.clip(0, a, b)
    endpoints = np.array([a, b])
    pa, pb = np.exp(-0.5 * (endpoints - anchor) * (endpoints + anchor) - mass)
    return float(sd**2 * (1 + a * pa - b * pb))


def fit_deviation(moment: float, lower: float, upper: float) -> float:
    """Match the truncated second moment, including its uniform limit."""
    if moment >= second_moment(np.inf, lower, upper):
        return np.inf
    if moment <= second_moment(MIN_DEVIATION, lower, upper):
        return MIN_DEVIATION
    lo, hi = np.log(MIN_DEVIATION), np.log(max(1.0, abs(lower), abs(upper)))
    for _ in range(128):
        if second_moment(np.exp(hi), lower, upper) >= moment:
            break
        hi += np.log(2)
    else:
        raise FloatingPointError("could not bracket truncated deviation")
    for _ in range(64):
        mid = (lo + hi) / 2
        if second_moment(np.exp(mid), lower, upper) < moment:
            lo = mid
        else:
            hi = mid
    return float(np.exp((lo + hi) / 2))


def component_objective(moment: float, sd: float, lower: float, upper: float) -> float:
    if np.isposinf(sd):
        return -np.log(upper - lower)
    anchor = np.clip(0, lower, upper)
    return float(
        -np.log(sd) - log_scaled_mass(lower / sd, upper / sd) - 0.5 * (moment - anchor**2) / sd / sd
    )


def maximization(
    values: FloatArray,
    responsibilities: FloatArray,
    means: FloatArray,
    deviations: FloatArray,
    bounds: FloatArray,
) -> tuple[FloatArray, FloatArray]:
    """Maximize the truncated objective with fixed means and sd >= 0.001."""
    if np.any(deviations < MIN_DEVIATION) or np.isnan(deviations).any():
        raise ValueError("fitted deviations must be at least 0.001")
    mass = np.sum(responsibilities, axis=0)
    if not np.isfinite(mass).all() or mass.sum() <= 0 or np.any(responsibilities < 0):
        raise ValueError("responsibilities must have finite positive total mass")
    new_deviations = deviations.copy()
    for j, (lower, upper) in enumerate(bounds):
        if not lower < upper:
            raise ValueError("bounds must have positive width")
        selected = responsibilities[:, j] > 0
        if not selected.any():
            continue
        x = values[selected]
        if np.any((x < lower) | (x > upper)):
            raise ValueError("responsibilities outside component support")
        moment = float(np.dot(responsibilities[selected, j] / mass[j], (x - means[j]) ** 2))
        lower, upper = lower - means[j], upper - means[j]
        candidate = fit_deviation(moment, lower, upper)
        old_q = component_objective(moment, deviations[j], lower, upper)
        new_q = component_objective(moment, candidate, lower, upper)
        if not np.isfinite([old_q, new_q]).all():
            raise FloatingPointError("non-finite truncated objective")
        if new_q >= old_q:
            new_deviations[j] = candidate
    return new_deviations, mass / mass.sum()
