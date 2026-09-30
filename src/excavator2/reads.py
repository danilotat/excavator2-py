"""MAPQ-filtered BAM selection and independent inclusive window counts."""

import numpy as np
import pysam

from .artifacts import coordinate_values

CHUNK_SIZE = 500000


def count_positions(chunks, starts, ends):
    starts, ends = np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64)
    if starts.ndim != 1 or starts.shape != ends.shape or not len(starts):
        raise ValueError("counting requires matching nonempty window vectors")
    if (starts > ends).any() or (np.diff(starts) < 0).any():
        raise ValueError("counting requires windows ordered by start")
    counts = np.zeros(len(starts), dtype=np.int64)
    previous = None
    for chunk in chunks:
        chunk = np.asarray(chunk, dtype=np.int64)
        if chunk.ndim != 1:
            raise ValueError("selected positions must be a vector")
        if not len(chunk):
            continue
        if (previous is not None and chunk[0] < previous) or (np.diff(chunk) < 0).any():
            raise ValueError("selected positions must be sorted")
        previous = chunk[-1]
        counts += np.searchsorted(chunk, ends, side="right") - np.searchsorted(
            chunk, starts, side="left"
        )
    return counts


def selected_chunks(bam, chromosome, mapq=20):
    if isinstance(mapq, bool) or not isinstance(mapq, (int, np.integer)) or not 0 <= mapq <= 255:
        raise ValueError("MAPQ must be an integer between 0 and 255")
    chunk = []
    for read in bam.fetch(chromosome):
        if read.flag & 1028 or read.mapping_quality < mapq:
            continue
        chunk.append(read.reference_start + 1)  # SAM POS, not pysam's zero-based start
        if len(chunk) == CHUNK_SIZE:
            yield np.array(chunk, dtype=np.int64)
            chunk = []
    yield np.array(chunk, dtype=np.int64)


def count_bam(path, target, chromosomes, mapq=20):
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
                    selected_chunks(bam, chromosome, mapq),
                    coordinate_values(windows[:, 1]),
                    coordinate_values(windows[:, 2]),
                )
            )
            selected_order.extend(indices)
    if not np.array_equal(selected_order, np.arange(len(target))):
        raise ValueError("target rows must follow the declared chromosome order")
    return np.concatenate(counts)
