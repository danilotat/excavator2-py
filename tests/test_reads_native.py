"""Native BAM counting agrees with independent per-window selection."""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pysam
import pytest
from numpy.testing import assert_array_equal

from excavator2 import _core
from excavator2.reads import count_bam


def native_count(path, chromosome, starts, ends, mapq):
    plan = _core.RegionPlan(chromosome, starts, ends)
    reader = _core.BamReader(path)
    try:
        return reader.count(plan, mapq)
    finally:
        reader.close()


def make_bam(path, positions, flags, qualities):
    with pysam.AlignmentFile(
        path, "wb", header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
    ) as bam:
        for i, (position, flag, quality) in enumerate(
            zip(positions, flags, qualities, strict=True)
        ):
            read = pysam.AlignedSegment()
            read.query_name = str(i)
            read.reference_id = 0
            read.reference_start = int(position) - 1
            read.flag = int(flag)
            read.mapping_quality = int(quality)
            read.query_sequence = "A"
            read.cigarstring = "1M"
            bam.write(read)
    pysam.index(str(path))


@pytest.mark.parametrize("mapq", [0, 20, 60, 255])
def test_native_overlaps_endpoints_flags_and_concurrent_calls(tmp_path, mapq):
    rng = np.random.default_rng(413)
    positions = np.sort(rng.integers(1, 1000, size=20000))
    flags = rng.choice([0, 1, 2, 4, 16, 64, 128, 256, 512, 1024, 2048], size=len(positions))
    qualities = rng.choice([0, 19, 20, 60, 255], size=len(positions))
    path = tmp_path / "reads.bam"
    make_bam(path, positions, flags, qualities)
    starts = np.sort(np.r_[0, 1, 1, rng.integers(1, 1100, size=200), 2000])
    ends = starts + rng.integers(0, 300, size=len(starts))
    selected = positions[((flags & 1028) == 0) & (qualities >= mapq)]
    expected = [np.count_nonzero((selected >= a) & (selected <= b)) for a, b in zip(starts, ends)]

    plan = _core.RegionPlan("chr1", starts, ends)

    def count(_):
        reader = _core.BamReader(str(path))
        try:
            return reader.count(plan, mapq)
        finally:
            reader.close()

    with ThreadPoolExecutor(max_workers=4) as pool:
        for actual in pool.map(count, range(8)):
            assert_array_equal(actual, expected)


def test_native_empty_and_terminal_windows(tmp_path):
    path = tmp_path / "empty.bam"
    make_bam(path, [], [], [])
    assert_array_equal(native_count(str(path), "chr1", [1], [20], 20), [0])
    path = tmp_path / "terminal.bam"
    make_bam(path, [10, 10, 20, 21], [0] * 4, [60] * 4)
    assert_array_equal(
        native_count(str(path), "chr1", [10, 10, 20, 21], [10, 20, 20, 21], 20),
        [2, 3, 1, 1],
    )
    with pytest.raises(ValueError, match="chromosome"):
        native_count(str(path), "missing", [1], [2], 20)


@pytest.mark.parametrize(
    "starts,ends,mapq", [([], [], 20), ([2, 1], [3, 3], 20), ([2], [1], 20), ([1], [2], 256)]
)
def test_native_rejects_invalid_arguments(starts, ends, mapq):
    with pytest.raises(ValueError):
        native_count("unused.bam", "chr1", starts, ends, mapq)


@pytest.mark.parametrize("mapq", [True, -1, 256, 20.5])
def test_python_boundary_rejects_invalid_mapq(mapq):
    with pytest.raises(ValueError, match="MAPQ"):
        count_bam("unused.bam", np.empty((0, 3)), [], mapq)


@pytest.mark.parametrize("threads", [1, 2, 4, 16])
def test_target_batches_preserve_overlaps_and_row_order(tmp_path, threads):
    path = tmp_path / "regions.bam"
    positions = np.arange(1, 2001)
    make_bam(path, positions, np.zeros(len(positions)), np.full(len(positions), 60))
    starts = np.arange(10, 1900, 10)
    # Nested, overlapping, and point windows cross task boundaries.
    ends = starts + np.where(np.arange(len(starts)) % 3 == 0, 100, 0)
    target = np.column_stack([np.repeat("chr1", len(starts)), starts, ends])
    expected = [np.count_nonzero((positions >= a) & (positions <= b)) for a, b in zip(starts, ends)]
    assert_array_equal(count_bam(path, target, ["chr1"], threads=threads), expected)


def test_index_query_does_not_count_reads_starting_before_window(tmp_path):
    path = tmp_path / "long-read.bam"
    with pysam.AlignmentFile(path, "wb", header={"SQ": [{"SN": "chr1", "LN": 1000}]}) as bam:
        read = pysam.AlignedSegment()
        read.query_name = "spans-window"
        read.reference_id = 0
        read.reference_start = 0
        read.mapping_quality = 60
        read.query_sequence = "A" * 200
        read.cigarstring = "200M"
        bam.write(read)
    pysam.index(str(path))
    assert_array_equal(native_count(str(path), "chr1", [100], [150], 20), [0])


def test_region_jobs_split_distant_windows():
    from excavator2.reads import region_jobs

    starts = np.r_[np.arange(100), np.arange(100) + 10_000_000]
    target = np.column_stack([np.repeat("chr1", len(starts)), starts, starts + 5])
    jobs = list(region_jobs(target, ["chr1"]))
    assert_array_equal(np.concatenate([job[0] for job in jobs]), np.arange(len(target)))
    assert all(starts[job[0][-1]] - starts[job[0][0]] <= 1_000_000 for job in jobs)


def test_reader_reuses_index_for_out_of_order_queries(tmp_path):
    path = tmp_path / "reads.bam"
    make_bam(path, [10, 20, 30], [0] * 3, [60] * 3)
    reader = _core.BamReader(str(path))
    plans = [_core.RegionPlan("chr1", [start], [start]) for start in [30, 10, 20]]
    for plan in plans * 3:
        assert_array_equal(reader.count(plan, 20), [1])
    with pytest.raises(ValueError, match="chromosome"):
        reader.count(_core.RegionPlan("missing", [10], [20]), 20)
    assert_array_equal(reader.count(plans[0], 20), [1])
    reader.close()
    reader.close()
    with pytest.raises(ValueError, match="closed"):
        reader.count(plans[0], 20)


def test_worker_reuses_and_releases_readers(tmp_path, monkeypatch):
    import weakref

    from excavator2.reads import RegionCounter

    paths = [tmp_path / f"sample-{i}.bam" for i in range(2)]
    for i, path in enumerate(paths):
        make_bam(path, [10] * (i + 1), [0] * (i + 1), [60] * (i + 1))
    factory = _core.BamReader
    opened, closed, references = [], [], []

    class TrackedReader:
        def __init__(self, path):
            self.reader = factory(path)
            self.path = path
            opened.append(path)
            references.append(weakref.ref(self))

        def count(self, plan, mapq):
            return self.reader.count(plan, mapq)

        def close(self):
            closed.append(self.path)
            self.reader.close()

    monkeypatch.setattr(_core, "BamReader", TrackedReader)
    counter = RegionCounter()
    job = (None, _core.RegionPlan("chr1", [10], [10]))
    with ThreadPoolExecutor(max_workers=1) as pool:
        for index in [0, 0, 1, 1]:
            assert_array_equal(pool.submit(counter, paths[index], job).result(), [index + 1])
        with pytest.raises(ValueError, match="open BAM"):
            pool.submit(counter, tmp_path / "missing.bam", job).result()
        assert_array_equal(pool.submit(counter, paths[0], job).result(), [1])
    assert opened == list(map(str, [paths[0], paths[1], paths[0]]))
    assert closed == list(map(str, paths))
    assert all(reference() is None for reference in references)
