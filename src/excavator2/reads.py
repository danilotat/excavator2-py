"""Parallel indexed BAM counting over batches of target windows."""

from concurrent.futures import ThreadPoolExecutor
from threading import local

import numpy as np

from . import _core
from .artifacts import coordinate_values


def region_jobs(target, chromosomes, threads=1):
    """Partition target rows into nearby batches, without splitting any window."""
    if threads < 1:
        raise ValueError("threads must be positive")
    if not chromosomes or len(set(chromosomes)) != len(chromosomes):
        raise ValueError("require nonempty unique target chromosomes")
    groups = [
        (chromosome, np.flatnonzero(target[:, 0] == chromosome)) for chromosome in chromosomes
    ]
    order = np.concatenate([indices for _, indices in groups])
    if not len(order) or not np.array_equal(order, np.arange(len(target))):
        raise ValueError("target rows must follow the declared chromosome order")
    batch_size = max(1, min(1024, (len(target) + 4 * threads - 1) // (4 * threads)))
    for chromosome, indices in groups:
        starts = coordinate_values(target[indices, 1])
        ends = coordinate_values(target[indices, 2])
        if not len(starts) or (starts > ends).any() or (np.diff(starts) < 0).any():
            raise ValueError("counting requires nonempty windows ordered by start")
        offset = 0
        while offset < len(indices):
            stop = min(offset + batch_size, len(indices))
            # Avoid scanning large gaps between small, distant target windows.
            nearby = int(np.searchsorted(starts, starts[offset] + 1_000_000, side="right"))
            stop = min(stop, nearby)
            yield (
                indices[offset:stop],
                _core.RegionPlan(str(chromosome), starts[offset:stop], ends[offset:stop]),
            )
            offset = stop


class RegionCounter:
    """Keep one BAM/index per worker, replacing it when the input changes."""

    def __init__(self):
        self._local = local()

    def __call__(self, path, job, mapq=20):
        if (
            isinstance(mapq, bool)
            or not isinstance(mapq, (int, np.integer))
            or not 0 <= mapq <= 255
        ):
            raise ValueError("MAPQ must be an integer between 0 and 255")
        path = str(path)
        cached = getattr(self._local, "cached", None)
        if cached is None or cached[0] != path:
            if cached is not None:
                cached[1].close()
                del self._local.cached
            cached = (path, _core.BamReader(path))
            self._local.cached = cached
        return cached[1].count(job[1], int(mapq))


def count_bam(path, target, chromosomes, mapq=20, threads=1):
    if isinstance(mapq, bool) or not isinstance(mapq, (int, np.integer)) or not 0 <= mapq <= 255:
        raise ValueError("MAPQ must be an integer between 0 and 255")
    jobs = list(region_jobs(target, chromosomes, threads))
    counts = np.zeros(len(target), dtype=np.int64)
    counter = RegionCounter()
    with ThreadPoolExecutor(max_workers=min(threads, len(jobs))) as pool:
        for job, result in zip(jobs, pool.map(lambda job: counter(path, job, mapq), jobs)):
            counts[job[0]] = result
    return counts
