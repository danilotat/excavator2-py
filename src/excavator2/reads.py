"""Legacy BAM selection and chunk-count policy, using NumPy's native searches.

The old awk predicate retains every MAPQ. Only flags 0x4 and 0x400 are excluded.
Counting follows MakeReadCount.R + EXOMECOUNT, including lost residuals and the
unflushed final interval. No custom C++ reader or counting kernel is needed.
"""

import numpy as np
import pysam

from .artifacts import coordinate_values

CHUNK_SIZE = 500000


def count_positions(chunks, starts, ends):
    starts, ends = np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64)
    if starts.ndim != 1 or starts.shape != ends.shape or not len(starts):
        raise ValueError("counting requires matching nonempty window vectors")
    if (starts > ends).any() or (np.diff(starts) < 0).any():
        raise ValueError("legacy counting requires windows ordered by start")
    effective_ends = np.maximum.accumulate(ends)
    counts = np.zeros(len(starts), dtype=np.int64)
    window, residual, previous = 0, 0, -1
    for chunk in chunks:
        chunk = np.asarray(chunk, dtype=np.int64)
        if chunk.ndim != 1 or not len(chunk):
            # The R wrapper fails reading an empty chunk, including exact multiples
            # of 500,000 unless it already advanced past the last target window.
            raise ValueError("legacy MakeReadCount cannot process an empty chunk")
        if chunk[0] < previous or (np.diff(chunk) < 0).any():
            raise ValueError("selected positions must be sorted")
        previous = chunk[-1]
        if chunk[-1] <= ends[window]:
            residual = len(chunk)  # assignment, not addition: preserve the wrapper
        else:
            stop = int(np.searchsorted(effective_ends, chunk[-1], side="left"))
            indices = np.arange(window, stop)
            upper = np.searchsorted(chunk, effective_ends[indices], side="right")
            lower = np.searchsorted(chunk, starts[indices], side="left")
            if len(indices) > 1:
                lower[1:] = np.maximum(lower[1:], upper[:-1])
            counts[indices] = upper - lower
            counts[window] += residual
            cursor = int(upper[-1])
            window = stop
            if window == len(starts):
                # The original dereferences beyond its Fortran arrays here. All
                # defined window counts have already been written; stop safely.
                break
            lower = max(cursor, int(np.searchsorted(chunk, starts[window], side="left")))
            residual = len(chunk) - lower
        if len(chunk) < CHUNK_SIZE:
            break
    return counts


def selected_chunks(bam, chromosome):
    chunk = []
    for read in bam.fetch(chromosome):
        if read.flag & 1028:
            continue
        chunk.append(read.reference_start + 1)  # SAM POS, not pysam's zero-based start
        if len(chunk) == CHUNK_SIZE:
            yield np.array(chunk, dtype=np.int64)
            chunk = []
    yield np.array(chunk, dtype=np.int64)


def count_bam(path, target, chromosomes):
    with pysam.AlignmentFile(str(path), "rb") as bam:
        if bam.is_cram:
            raise ValueError("CRAM preparation is not yet qualified")
        if not bam.has_index():
            raise ValueError("legacy-compatible BAM preparation requires an index")
        counts = []
        selected_order = []
        for chromosome in chromosomes:
            indices = np.flatnonzero(target[:, 0] == chromosome)
            windows = target[indices]
            counts.append(
                count_positions(
                    selected_chunks(bam, chromosome),
                    coordinate_values(windows[:, 1]),
                    coordinate_values(windows[:, 2]),
                )
            )
            selected_order.extend(indices)
    if not np.array_equal(selected_order, np.arange(len(target))):
        raise ValueError("target rows must follow the declared chromosome order")
    return np.concatenate(counts)
