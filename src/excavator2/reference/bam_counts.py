"""Python counting oracle for native BAM counting comparisons."""

import numpy as np

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
