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
    # Contigs outside autosomes and X are reported; centromere names may omit "chr".
    with (tmp_path / "target.bed").open("a") as handle:
        handle.write("chrY\t800\t850\nchrM\t10\t60\nchrY\t900\t950\n")
    centro = tmp_path / "centro.tsv"
    centro.write_text(centro.read_text().replace("chr", ""))
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
    assert "excavator2 target: [1/4] Loading" in result.stderr
    assert "excavator2 target: [3/4] Extracting reference features" in result.stderr
    assert "excavator2 target: Complete:" in result.stderr
    assert "outside reference autosomes and X: chrM (1), chrY (2)" in result.stderr
    manifest = load_manifest(output, "target")
    assert manifest["excluded_bed_contigs"] == {"chrM": 1, "chrY": 2}
    assert manifest["flank"] == 200
    assert manifest["chromosomes"] == [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    assert manifest["centromeres"]["chrX"] == [1400, 1600]
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


def test_target_flank_setting_moves_off_target_windows(config, tmp_path):
    settings = yaml.safe_load(config.read_text())
    settings["Target"]["Flank"] = 300
    config.write_text(yaml.safe_dump(settings))
    run_target(config, tmp_path / "target")
    manifest = load_manifest(tmp_path / "target", "target")
    assert manifest["flank"] == 300
    with np.load(tmp_path / "target" / manifest["preparation"], allow_pickle=False) as data:
        target = data["target"]
    outside = target[(target[:, 0] == "chr1") & (target[:, 4] == "OUT"), 1:3].astype(int)
    # Flanks leave no window before the first target (500-560).
    assert outside[0, 0] == 861
    assert np.all((outside[:, 1] <= 1000) | (outside[:, 0] > 1700))


def test_target_failure_does_not_publish(config, tmp_path):
    settings = yaml.safe_load(config.read_text())
    Path(settings["Reference"]["Gaps"]).write_text("id\tchrom\tstart\tend\n0\tchr1\t0\t1\n")
    with pytest.raises(ValueError, match="alternate-contig gap row"):
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


@pytest.mark.parametrize("failure", ["features", "write"])
def test_target_reference_and_write_errors_are_atomic(config, tmp_path, monkeypatch, failure):
    if failure == "features":
        settings = yaml.safe_load(config.read_text())
        coordinates = Path(settings["Reference"]["Chromosomes"])
        coordinates.write_text(coordinates.read_text().replace("3000", "6000"))
        message = "invalid reference window"
    else:

        def fail_write(*args, **kwargs):
            raise OSError("injected disk failure")

        monkeypatch.setattr(np, "savez_compressed", fail_write)
        message = "injected disk failure"
    with pytest.raises((ValueError, OSError), match=message):
        run_target(config, tmp_path / "failed")
    assert not (tmp_path / "failed").exists()
    assert not list(tmp_path.glob(".excavator2-target-*"))
