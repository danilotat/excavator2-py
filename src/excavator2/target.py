"""Target geometry with named bounds, merged targets and bounded flanks."""

from pathlib import Path

import numpy as np


def target_geometry(bed: Path, coordinates: Path, gaps: Path, window: int) -> np.ndarray:
    """Return chromosome, start, end, ID and IN/OUT columns."""
    if isinstance(window, bool) or not isinstance(window, (int, np.integer)) or window < 10:
        raise ValueError("window must be an integer of at least 10")
    bed_rows = np.loadtxt(bed, dtype=str, delimiter="\t", usecols=(0, 1, 2), ndmin=2)
    coord_rows = np.loadtxt(coordinates, dtype=str, delimiter="\t", usecols=(0, 1, 2), ndmin=2)
    gap_rows = np.loadtxt(gaps, dtype=str, delimiter="\t", skiprows=1, usecols=(1, 2, 3), ndmin=2)
    if len(bed_rows) == 0 or len(coord_rows) < 23:
        raise ValueError("geometry requires a nonempty BED and 23 canonical coordinate rows")
    bed_positions = bed_rows[:, 1:3].astype(np.int64)
    bounds = coord_rows[:, 1:3].astype(np.int64)
    gap_positions = gap_rows[:, 1:3].astype(np.int64)
    if np.any(bounds[:, 1] <= bounds[:, 0]):
        raise ValueError("chromosome bounds must have positive length")
    if np.any(gap_positions[:, 1] < gap_positions[:, 0]):
        raise ValueError("gap ends must not precede starts")

    prefix = "chr" if len(bed_rows[0, 0]) > 3 else ""
    chromosomes = [f"{prefix}{c}" for c in [*range(1, 23), "X"]]
    gap_names = np.array([c if prefix else c.replace("chr", "", 1) for c in gap_rows[:, 0]])
    alternate = np.array(["_" in c for c in gap_names], dtype=bool)
    # R's -integer(0) selects zero rows: preserve the failure, not a silent fix.
    if not np.any(alternate):
        raise ValueError("legacy gap filtering requires an alternate-contig gap row")
    gap_names, gap_positions = gap_names[~alternate], gap_positions[~alternate]

    named_bounds = {}
    for name, bound in zip(coord_rows[:, 0], bounds, strict=True):
        name = name.removeprefix("chr")
        if name in named_bounds:
            raise ValueError(f"duplicate chromosome bounds for {name}")
        named_bounds[name] = bound
    if any(c.removeprefix("chr") not in named_bounds for c in chromosomes):
        raise ValueError("geometry requires bounds for all 23 canonical chromosomes")

    result = []
    for chromosome in chromosomes:
        start, end = named_bounds[chromosome.removeprefix("chr")]
        inside = bed_positions[bed_rows[:, 0] == chromosome]
        if len(inside) and (
            np.any(inside[:, 0] < start)
            or np.any(inside[:, 1] > end)
            or np.any(inside[:, 1] <= inside[:, 0])
        ):
            raise ValueError(f"invalid target interval for {chromosome}")
        merged = []
        for left, right in sorted(map(tuple, inside)):
            if merged and left <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], right)
            else:
                merged.append([left, right])
        inside = np.asarray(merged, dtype=np.int64).reshape(-1, 2)
        outside = []
        if len(inside):
            empty_starts = np.concatenate(([start], inside[:, 1]))
            empty_ends = np.concatenate((inside[:, 0], [end]))
            for left, right in zip(empty_starts, empty_ends, strict=True):
                left, right = left + 200, right - 200
                boundaries = left + np.arange(max(0, (right - left) // window) + 1) * window
                outside.extend(zip(boundaries[:-1] + 1, boundaries[1:], strict=True))
        else:
            boundaries = np.arange(start, end + 1, window)
            outside.extend(zip(boundaries[:-1] + 1, boundaries[1:], strict=True))
        outside = np.asarray(outside, dtype=np.int64).reshape(-1, 2)
        positions = np.concatenate((inside, outside))
        labels = np.array(["IN"] * len(inside) + ["OUT"] * len(outside))
        order = np.argsort(positions[:, 0], kind="stable")
        positions, labels = positions[order], labels[order]
        chromosome_gaps = gap_positions[gap_names == chromosome]
        if not len(chromosome_gaps):
            raise ValueError(f"legacy gap filtering requires gaps for {chromosome}")
        keep = np.ones(len(positions), dtype=bool)
        # Vectorize across windows without allocating a windows-by-gaps matrix.
        for left, right in chromosome_gaps:
            keep &= ~np.any((positions >= left) & (positions < right), axis=1)
        for (left, right), label in zip(positions[keep], labels[keep], strict=True):
            result.append((chromosome, str(left), str(right), f"a{len(result) + 1}", label))
    if not result:
        raise ValueError("legacy geometry has no surviving windows")
    return np.asarray(result, dtype=str)
