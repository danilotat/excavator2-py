"""Literal size → MAP → GC median corrections from LibraryExomeRC.R."""

import numpy as np

from .artifacts import coordinate_values


def correct(values, covariate, maximum=100):
    values, covariate = np.asarray(values, dtype=float), np.asarray(covariate, dtype=float)
    if values.ndim != 1 or values.shape != covariate.shape or not len(values):
        raise ValueError("normalization requires matching nonempty vectors")
    if not np.isfinite(values).all() or not np.isfinite(covariate).all():
        raise ValueError("non-finite normalization input")
    result = values.copy()
    master = np.median(values)
    # R seq(0, max(length), by=5) stops before a partial final size bin.
    for upper in range(5, int(maximum) + 1, 5):
        selected = (covariate <= upper) & (covariate >= 0 if upper == 5 else covariate > upper - 5)
        if selected.any():
            median = np.median(values[selected])
            if median > 0:
                result[selected] = values[selected] * master / median
    return result


def normalize(counts, target, gc, mappability):
    counts, gc, mappability = (np.asarray(x, dtype=float) for x in [counts, gc, mappability])
    n = len(target)
    if any(x.shape != (n,) for x in [counts, gc, mappability]):
        raise ValueError("counts and target covariates have different shapes")
    length = coordinate_values(target[:, 2]) - coordinate_values(target[:, 1])
    if (length <= 0).any() or (counts < 0).any():
        raise ValueError("target lengths must be positive and counts nonnegative")
    if not np.isin(target[:, 4], ["IN", "OUT"]).all():
        raise ValueError("unknown target class")
    weighted = counts / length
    size, mapped, corrected, final = (weighted.copy() for _ in range(4))
    for kind in ["IN", "OUT"]:
        selected = target[:, 4] == kind
        if not selected.any():
            raise ValueError("legacy normalization requires IN and OUT windows")
        if kind == "IN":
            if max(length[selected]) < 5:
                raise ValueError("legacy size correction requires at least one complete bin")
            size[selected] = correct(weighted[selected], length[selected], max(length[selected]))
        mapped[selected] = correct(size[selected], mappability[selected] * 100)
        corrected[selected] = correct(mapped[selected], gc[selected] * 100)
        values = corrected[selected].copy()
        nonzero = values[values != 0]
        if not len(nonzero):
            raise ValueError(f"legacy normalization has no nonzero {kind} counts")
        values[values == 0] = min(nonzero)
        final[selected] = values
    if not np.isfinite(final).all():
        raise ValueError("legacy normalization produced non-finite counts")
    return {"weighted": weighted, "size": size, "map": mapped, "gc": corrected, "normalized": final}
