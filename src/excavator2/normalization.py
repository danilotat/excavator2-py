"""Preserve observed depth; calibration belongs to test/reference analysis."""

import numpy as np

from .artifacts import coordinate_values

DEPTH_POLICY = "raw-counts-v1"


def normalize(counts, target, gc, mappability):
    """Return uncorrected read density without imputing zero observations."""
    counts, gc, mappability = (np.asarray(x, dtype=float) for x in [counts, gc, mappability])
    n = len(target)
    if not n or any(x.shape != (n,) for x in [counts, gc, mappability]):
        raise ValueError("counts and target covariates require matching nonempty shapes")
    if any(not np.isfinite(x).all() for x in [counts, gc, mappability]):
        raise ValueError("missing or non-finite preparation input")
    length = coordinate_values(target[:, 2]) - coordinate_values(target[:, 1])
    if (length <= 0).any() or (counts < 0).any() or (counts != np.floor(counts)).any():
        raise ValueError("target lengths must be positive and counts nonnegative integers")
    if not np.isin(target[:, 4], ["IN", "OUT"]).all():
        raise ValueError("unknown target class")
    weighted = counts / length
    return {"weighted": weighted, "normalized": weighted.copy()}
