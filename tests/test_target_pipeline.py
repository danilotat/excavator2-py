"""Target CLI publication, artifact compatibility and failure handling."""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pyBigWig
import pysam
import pytest
import yaml
from numpy.testing import assert_array_equal

from excavator2.artifacts import load_manifest, metadata_identity, window_metadata
from excavator2.target_pipeline import run_target


@pytest.fixture
def config(tmp_path):
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    fasta = tmp_path / "reference.fa"
    fasta.write_text("".join(f">{c}\n{'G' * 5000}\n" for c in chromosomes))
    pysam.faidx(str(fasta))
    bigwig = tmp_path / "reference.bw"
    with pyBigWig.open(str(bigwig), "w") as bw:
        bw.addHeader([(c, 5000) for c in chromosomes])
        bw.addEntries(chromosomes, [0] * 23, ends=[5000] * 23, values=[0.2] * 23)
    coords, gaps, centro, bed = [
        tmp_path / name for name in ["coords.tsv", "gaps.tsv", "centro.tsv", "target.bed"]
    ]
    coords.write_text("".join(f"{c}\t0\t3000\n" for c in chromosomes))
    # Remove each zero-start window; the original extraction rejects negative BED starts.
    gaps.write_text(
        "id\tchrom\tstart\tend\n"
        + "".join(f"0\t{c}\t0\t1\n" for c in chromosomes)
        + "0\tchr1_alt\t0\t1\n"
    )
    centro.write_text("chrom\tstart\tend\n" + "".join(f"{c}\t1400\t1600\n" for c in chromosomes))
    bed.write_text("chr1\t500\t560\nchr1\t1300\t1400\n")
    settings = {
        "Reference": {
            "Assembly": "synthetic",
            "FASTA": str(fasta),
            "BigWig": str(bigwig),
            "Chromosomes": str(coords),
            "Gaps": str(gaps),
            "Centromeres": str(centro),
        },
        "Target": {"Name": "panel", "BED": str(bed), "Window": 100},
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(settings))
    return path


def test_target_cli_publishes_compatible_artifact(config, tmp_path):
    output = tmp_path / "target"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "excavator2",
            "target",
            "--config",
            str(config),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    manifest = load_manifest(output, "target")
    with np.load(output / manifest["preparation"], allow_pickle=False) as data:
        target = data["target"]
        assert manifest["target_id"] == metadata_identity(window_metadata(target))
        assert_array_equal(data["gc"], np.ones(len(target)))
        assert_array_equal(data["mappability"], np.full(len(target), 0.2))
    assert_array_equal(np.loadtxt(output / "Filtered.txt", dtype=str), target)
    for chrom in manifest["chromosomes"]:
        with np.load(output / manifest["references"][chrom]) as data:
            assert_array_equal(data["matrix"][:, 0], target[target[:, 0] == chrom, 1])
            assert set(data["matrix"][:, 1]) == {"G"}
    assert len(manifest["inputs"]) == 6
    assert not list(tmp_path.glob(".excavator2-target-*"))


def test_target_failure_does_not_publish(config, tmp_path):
    settings = yaml.safe_load(config.read_text())
    Path(settings["Reference"]["Gaps"]).write_text("id\tchrom\tstart\tend\n0\tchr1_alt\t0\t1\n")
    with pytest.raises(ValueError, match="requires gaps"):
        run_target(config, tmp_path / "failed")
    assert not (tmp_path / "failed").exists()
    assert not list(tmp_path.glob(".excavator2-target-*"))


def test_target_refuses_existing_output_and_force(config, tmp_path):
    output = tmp_path / "target"
    output.mkdir()
    sentinel = output / "keep.txt"
    sentinel.write_text("keep")
    with pytest.raises(ValueError, match="new output directory"):
        run_target(config, output)
    with pytest.raises(ValueError, match="force is not supported"):
        run_target(config, tmp_path / "forced", force=True)
    assert sentinel.read_text() == "keep"
    assert not (tmp_path / "forced").exists()
