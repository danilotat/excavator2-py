"""Legacy target geometry, separate from reference features and artifact writing.

Coordinates deliberately retain FilterTarget.R's conventions and known quirks.
See docs/target-porting.md before changing ordering, flanks or gap filtering.
"""

from pathlib import Path

import numpy as np


def target_geometry(bed: Path, coordinates: Path, gaps: Path, window: int) -> np.ndarray:
    """Return the legacy five-column string matrix (chrom, start, end, ID, IN/OUT).

    Input files follow the original tab-separated BED, chromosome coordinates and
    header-bearing gap formats. Only integer genomic coordinates are supported.
    Reference feature extraction and output-directory creation are not performed.
    """
    if isinstance(window, bool) or not isinstance(window, (int, np.integer)) or window < 10:
        raise ValueError("window must be an integer of at least 10")
    bed_rows = np.loadtxt(bed, dtype=str, delimiter="\t", usecols=(0, 1, 2), ndmin=2)
    coord_rows = np.loadtxt(coordinates, dtype=str, delimiter="\t", usecols=(1, 2), ndmin=2)
    gap_rows = np.loadtxt(gaps, dtype=str, delimiter="\t", skiprows=1, usecols=(1, 2, 3), ndmin=2)
    if len(bed_rows) == 0 or len(coord_rows) < 23:
        raise ValueError("geometry requires a nonempty BED and 23 canonical coordinate rows")
    bed_positions = bed_rows[:, 1:3].astype(np.int64)
    bounds = coord_rows.astype(np.int64)
    gap_positions = gap_rows[:, 1:3].astype(np.int64)
    if np.any(bounds[:23, 1] <= bounds[:23, 0]):
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

    result = []
    for index, chromosome in enumerate(chromosomes):
        start, end = bounds[index]
        inside = bed_positions[bed_rows[:, 0] == chromosome]
        outside = []
        if len(inside):
            # Do not sort or merge BED rows before constructing these gaps.
            empty_starts = np.concatenate(([start], inside[:, 1]))
            empty_ends = np.concatenate((inside[:, 0], [end]))
            eligible = empty_ends - empty_starts >= window + 400
            if not np.any(eligible):
                raise ValueError(f"legacy geometry has no eligible off-target gap for {chromosome}")
            for left, right in zip(empty_starts[eligible], empty_ends[eligible], strict=True):
                boundaries = left + 200 + np.arange((right - left) // window + 1) * window
                outside.append(np.column_stack((boundaries[:-1] + 1, boundaries[1:])))
        else:
            boundaries = np.arange(start, end + 1, window)
            outside.append(np.column_stack((boundaries[:-1], boundaries[1:] - 1)))
        outside = np.concatenate(outside)
        positions = np.concatenate((inside, outside))
        if not len(positions):
            raise ValueError(f"legacy geometry has no windows for {chromosome}")
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
