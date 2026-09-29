"""Exact counting and continuous normalization comparisons against original code."""

from pathlib import Path

import numpy as np
import pysam
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2.normalization import normalize
from excavator2.reads import CHUNK_SIZE, count_positions, selected_chunks

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
def test_counts_match_original_chunk_policy(name):
    with np.load(FIXTURES / f"{name}.npz") as f:
        reads = f["reads"]
        chunks = [reads[i : i + CHUNK_SIZE] for i in range(0, len(reads), CHUNK_SIZE)]
        if len(reads) % CHUNK_SIZE == 0:
            chunks.append(np.array([], dtype=int))
        if f["failed"][0]:
            with pytest.raises(ValueError, match="empty chunk"):
                count_positions(chunks, f["starts"], f["ends"])
        else:
            assert_array_equal(count_positions(chunks, f["starts"], f["ends"]), f["counts"])


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


def test_selection_keeps_mapq_secondary_supplementary_and_mates(tmp_path):
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
            if not flag & 1028:
                expected.append(i * 10 + 1)
    pysam.index(str(path))
    with pysam.AlignmentFile(path) as bam:
        assert_array_equal(np.concatenate(list(selected_chunks(bam, "chr1"))), expected)


def test_empty_stream_fails_and_final_window_is_flushed_only_by_later_read():
    with pytest.raises(ValueError, match="empty chunk"):
        count_positions([np.array([], dtype=int)], [10], [20])
    assert_array_equal(count_positions([[15, 20]], [10], [20]), [0])
    assert_array_equal(count_positions([[15, 20, 21]], [10], [20]), [2])


def test_numpy_counter_matches_literal_state_machine_on_overlapping_windows():
    from excavator2.reference.reads import count_positions as reference

    rng = np.random.default_rng(413)
    for _ in range(30):
        starts = np.sort(rng.integers(1, 1000, size=50))
        ends = starts + rng.integers(1, 100, size=50)
        reads = np.sort(rng.integers(1, 1500, size=1000))
        assert_array_equal(count_positions([reads], starts, ends), reference([reads], starts, ends))


@pytest.mark.parametrize("threads", [1, 4])
def test_prepare_cli_matches_legacy_and_is_thread_independent(tmp_path, threads):
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
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["threads"] == threads
    with np.load(fixture / "expected.npz") as old:
        for filename in manifest["samples"].values():
            with np.load(output / filename) as new:
                assert_array_equal(new["counts"], old["counts"])
                assert_array_equal(new["matrix"], old["matrix"])
                for field, legacy in [
                    ("weighted", "weighted"),
                    ("size", "size"),
                    ("map", "mapped"),
                    ("gc", "corrected"),
                    ("normalized", "normalized"),
                ]:
                    assert_allclose(new[field], old[legacy], rtol=1e-13, atol=2e-15)


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
