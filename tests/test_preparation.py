"""Exact counting and continuous normalization comparisons against original code."""

from pathlib import Path

import numpy as np
import pysam
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2.normalization import normalize
from excavator2.reference.bam_counts import CHUNK_SIZE, count_positions, selected_chunks

FIXTURES = Path(__file__).parent / "fixtures/legacy-preparation"


@pytest.mark.parametrize(
    "name",
    [
        "endpoints",
        "overlap",
        "final_only",
        "chunk_499999",
        "chunk_500000",
        "chunk_500001",
        "chunk_1000001",
    ],
)
def test_counts_correct_original_chunk_defects(name):
    with np.load(FIXTURES / f"{name}.npz") as f:
        reads = f["reads"]
        chunks = [reads[i : i + CHUNK_SIZE] for i in range(0, len(reads), CHUNK_SIZE)]
        if len(reads) % CHUNK_SIZE == 0:
            chunks.append(np.array([], dtype=int))
        expected = [
            np.count_nonzero((reads >= a) & (reads <= b))
            for a, b in zip(f["starts"], f["ends"], strict=True)
        ]
        assert_array_equal(count_positions(chunks, f["starts"], f["ends"]), expected)


def test_normalization_trace_matches_original():
    with np.load(FIXTURES / "normalization.npz") as f:
        n = len(f["counts"])
        target = np.column_stack(
            [
                np.repeat("chr1", n),
                np.zeros(n, dtype=int),
                f["lengths"],
                np.repeat("gene", n),
                f["classes"],
            ]
        )
        result = normalize(f["counts"], target, f["gc"], f["mappability"])
        for name, original in [
            ("weighted", "weighted"),
            ("size", "size"),
            ("map", "mapped"),
            ("gc", "corrected"),
            ("normalized", "normalized"),
        ]:
            assert_allclose(result[name], f[original], rtol=1e-14, atol=1e-15)


def test_selection_filters_mapq_and_retains_other_flag_semantics(tmp_path):
    path = tmp_path / "flags.bam"
    flags = [0, 1, 2, 4, 16, 64, 128, 256, 512, 1024, 2048]
    expected = []
    with pysam.AlignmentFile(
        path, "wb", header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    ) as bam:
        for i, flag in enumerate(flags):
            read = pysam.AlignedSegment()
            read.query_name = f"read{i}"
            read.flag = flag
            read.reference_id = 0
            read.reference_start = i * 10
            read.mapping_quality = [0, 1, 19, 20, 60][i % 5]
            read.query_sequence = "A"
            read.cigarstring = "1M"
            bam.write(read)
            if not flag & 1028 and read.mapping_quality >= 20:
                expected.append(i * 10 + 1)
    pysam.index(str(path))
    with pysam.AlignmentFile(path) as bam:
        assert_array_equal(np.concatenate(list(selected_chunks(bam, "chr1"))), expected)


def test_empty_stream_and_terminal_window():
    assert_array_equal(count_positions([[]], [10], [20]), [0])
    assert_array_equal(count_positions([[15, 20]], [10], [20]), [2])
    assert_array_equal(count_positions([[15, 20, 21]], [10], [20]), [2])


def test_counts_are_independent_of_partition_and_overlap():
    rng = np.random.default_rng(413)
    for _ in range(30):
        starts = np.sort(rng.integers(1, 1000, size=50))
        ends = starts + rng.integers(1, 100, size=50)
        reads = np.sort(rng.integers(1, 1500, size=1000))
        expected = [
            np.count_nonzero((reads >= a) & (reads <= b)) for a, b in zip(starts, ends, strict=True)
        ]
        for chunks in ([reads], np.array_split(reads, 17), [[], reads, []]):
            assert_array_equal(count_positions(chunks, starts, ends), expected)


@pytest.mark.parametrize("threads", [1, 4])
def test_prepare_cli_preserves_counts_and_is_thread_independent(tmp_path, threads):
    import json
    import subprocess
    import sys

    fixture = FIXTURES / "pipeline"
    samples = tmp_path / "samples.yaml"
    samples.write_text(f"Test1: {fixture / 'Test1.bam'}\nReplica: {fixture / 'Test1.bam'}\n")
    output = tmp_path / "prepared"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "excavator2",
            "prepare",
            "--samples",
            str(samples),
            "--target",
            str(fixture / "target"),
            "--output",
            str(output),
            "--threads",
            str(threads),
            "--mapq",
            "60",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "excavator2 prepare: Preparing 2 samples" in result.stderr
    assert "1/2 Prepared" in result.stderr
    assert "2/2 Prepared" in result.stderr
    assert "excavator2 prepare: Complete:" in result.stderr
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["threads"] == threads
    with np.load(fixture / "expected.npz") as old:
        for filename in manifest["samples"].values():
            with np.load(output / filename) as new:
                with pysam.AlignmentFile(fixture / "Test1.bam") as bam:
                    expected = []
                    for chrom, start, end in old["matrix"][:, [0, 2, 3]]:
                        expected.append(
                            sum(
                                not r.flag & 1028
                                and r.mapping_quality >= 60
                                and int(start) <= r.reference_start + 1 <= int(end)
                                for r in bam.fetch(chrom)
                            )
                        )
                assert_array_equal(new["counts"], expected)
                assert manifest["depth_policy"] == "size-map-gc-normalized-v1"
                assert_array_equal(
                    new["matrix"][:, [0, 1, 2, 3, 4, 6]], old["matrix"][:, [0, 1, 2, 3, 4, 6]]
                )
                assert_allclose(
                    new["weighted"],
                    np.array(expected)
                    / (old["matrix"][:, 3].astype(float) - old["matrix"][:, 2].astype(float)),
                )
                assert_allclose(new["matrix"][:, 5].astype(float), new["normalized"], atol=1e-7)


def test_failed_preparation_does_not_publish_output(tmp_path):
    import shutil

    from excavator2.prepare import run_preparation

    fixture = FIXTURES / "pipeline"
    bam = tmp_path / "no-index.bam"
    shutil.copyfile(fixture / "Test1.bam", bam)
    samples = tmp_path / "samples.yaml"
    samples.write_text(f"Test1: {bam}\n")
    with pytest.raises(ValueError, match="index"):
        run_preparation(samples, fixture / "target", tmp_path / "prepared")
    assert not (tmp_path / "prepared").exists()
    assert not list(tmp_path.glob(".excavator2-prepare-*"))


@pytest.mark.parametrize("threads", [1, 2])
def test_cli_worker_failure_is_nonzero_and_atomic(tmp_path, threads):
    import shutil
    import subprocess
    import sys

    fixture = FIXTURES / "pipeline"
    bad = tmp_path / "no-index.bam"
    shutil.copyfile(fixture / "Test1.bam", bad)
    samples = tmp_path / "samples.yaml"
    samples.write_text(f"Good: {fixture / 'Test1.bam'}\nBroken: {bad}\n")
    output = tmp_path / "prepared"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "excavator2",
            "prepare",
            "--samples",
            str(samples),
            "--target",
            str(fixture / "target"),
            "--output",
            str(output),
            "--threads",
            str(threads),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "Broken" in result.stderr and "index" in result.stderr
    assert "Complete:" not in result.stderr
    assert not output.exists()
    assert not list(tmp_path.glob(".excavator2-prepare-*"))


def test_one_sample_uses_region_workers_without_input_hashes(tmp_path, monkeypatch):
    import json
    import threading

    from excavator2 import prepare

    fixture = FIXTURES / "pipeline"
    samples = tmp_path / "samples.yaml"
    samples.write_text(f"Test1: {fixture / 'Test1.bam'}\n")
    output = tmp_path / "prepared"
    barrier = threading.Barrier(4)
    lock = threading.Lock()
    entered = 0
    thread_ids = set()
    original_count = prepare.RegionCounter.__call__
    original_digest = prepare.digest

    def count(*args):
        nonlocal entered
        with lock:
            entered += 1
            first_wave = entered <= 4
            thread_ids.add(threading.get_ident())
        if first_wave:
            barrier.wait(timeout=10)
        return original_count(*args)

    def digest(path):
        assert Path(path).suffix == ".npz", "preparation must not hash input files"
        return original_digest(path)

    monkeypatch.setattr(prepare.RegionCounter, "__call__", count)
    monkeypatch.setattr(prepare, "digest", digest)
    prepare.run_preparation(samples, fixture / "target", output, threads=4)
    assert len(thread_ids) == 4
    manifest = json.loads((output / "manifest.json").read_text())
    assert not {"bam_sha256", "target_manifest_sha256", "samples_sha256"} & manifest.keys()
    serial = tmp_path / "serial"
    monkeypatch.setattr(prepare.RegionCounter, "__call__", original_count)
    prepare.run_preparation(samples, fixture / "target", serial, threads=1)
    with np.load(output / "sample-0.npz") as parallel, np.load(serial / "sample-0.npz") as single:
        for name in parallel.files:
            assert_array_equal(parallel[name], single[name])


def test_boundary_plans_are_built_once_for_all_samples(tmp_path, monkeypatch):
    from excavator2 import _core, prepare

    fixture = FIXTURES / "pipeline"
    samples = tmp_path / "samples.yaml"
    samples.write_text(f"First: {fixture / 'Test1.bam'}\nSecond: {fixture / 'Test1.bam'}\n")
    factory = _core.RegionPlan
    plans = []
    uses = []
    count = prepare.RegionCounter.__call__

    def plan(*args):
        result = factory(*args)
        plans.append(result)
        return result

    def tracked_count(self, path, job, mapq):
        uses.append(job[1])
        return count(self, path, job, mapq)

    monkeypatch.setattr(_core, "RegionPlan", plan)
    monkeypatch.setattr(prepare.RegionCounter, "__call__", tracked_count)
    prepare.run_preparation(samples, fixture / "target", tmp_path / "prepared", threads=4)
    assert plans
    assert len(uses) == 2 * len(plans)
    assert all(sum(used is plan for used in uses) == 2 for plan in plans)
