"""Exercise the installed CLI against non-normal original-R paired/pooling outputs."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2.analyze import experimental_design, run_analysis, summarize
from excavator2.artifacts import digest, load_manifest, read_yaml

FIXTURE = Path(__file__).parent / "fixtures" / "legacy-analysis"


@pytest.fixture
def prepared_data():
    return FIXTURE / "inputs/prepared"


@pytest.mark.parametrize("experiment", ["paired", "pooling"])
def test_cli_reports_segment_interval_probabilities(tmp_path, experiment, prepared_data):
    output = tmp_path / "results"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "excavator2",
            "analyze",
            "--samples",
            str(FIXTURE / "samples.yaml"),
            "--input",
            str(prepared_data),
            "--target",
            str(FIXTURE / "inputs/target"),
            "--output",
            str(output),
            "--experiment",
            experiment,
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["fastcall_posterior"] == "truncated"
    for sample in ("Test1", "Test2"):
        folder = output / "Results" / sample
        hslm = f"HSLMResults_{sample}.txt"
        assert manifest["centering"] == "separate-IN-OUT-medians"
        assert "calibration" not in manifest
        assert "calibration_sha256" not in manifest
        assert "pseudocount_reads" not in manifest
        with np.load(folder / "checkpoints.npz") as data:
            expected = np.loadtxt(
                FIXTURE / "expected" / experiment / sample / hslm, dtype=str, skiprows=1
            )
            assert_allclose(data["ratios"], expected[:, 4].astype(float), atol=2e-8)
            table = np.loadtxt(folder / hslm, dtype=str, skiprows=1)
            values = table[data["segment_starts"], 5].astype(float)
            labels = np.searchsorted([-1.3, -0.5, 0.35, 0.9], values) - 2
            assert_array_equal(data["labels"], labels)
            assert_array_equal(data["posterior"], np.eye(5)[labels + 2])
            probabilities = data["call_probabilities"]
            assert_allclose(probabilities.sum(axis=1), 1)
            significant = labels != 0
            called = probabilities[np.arange(len(labels)), labels + 2][significant]
            altered = 1 - probabilities[significant, 2]
        assert manifest["fastcall_probability_semantics"] == "segment-level-interval-probability"
        assert set(manifest["samples"][sample]["measurement_noise"]) == {"IN", "OUT"}
        calls = np.genfromtxt(
            folder / f"FastCallResults_{sample}.txt", names=True, dtype=None, encoding="utf-8"
        )
        assert_allclose(calls["ProbCall"], called, atol=1e-14)
        assert_allclose(calls["ProbAltered"], altered, atol=1e-14)
        # Segments just below the gain/amplification edge: either class, surely altered.
        edge = np.isclose(calls["CorrectedSegment"], 0.9)
        assert edge.any()
        assert_allclose(calls["ProbCall"][edge], 0.5, atol=1e-6)
        assert_allclose(calls["ProbAltered"][edge], 1)


def test_failed_analysis_does_not_publish_partial_results(tmp_path, prepared_data):
    parameters = tmp_path / "parameters.yaml"
    parameters.write_text(
        "HSLM: {Omega: 0, Theta: 0.1, D_norm: 100}\n"
        "FastCall: {Cellularity: 1, d: 0.5, u: 0.35, minExons: 4}\n"
    )
    output = tmp_path / "results"
    with pytest.raises(ValueError, match="omega"):
        run_analysis(
            FIXTURE / "samples.yaml",
            prepared_data,
            FIXTURE / "inputs/target",
            output,
            "paired",
            parameters,
        )
    assert not output.exists()
    assert not list(tmp_path.glob(".excavator2-*"))


def test_truncated_support_failure_does_not_publish_partial_samples(
    tmp_path, monkeypatch, prepared_data
):
    from excavator2 import analyze

    original_segment_profile = analyze.segment_profile
    processed = []

    def unsupported_second_sample(*args):
        rows, path, indices, ids = original_segment_profile(*args)
        processed.append(True)
        if len(processed) == 2:
            rows[0, 5] = "100"
        return rows, path, indices, ids

    monkeypatch.setattr(analyze, "segment_profile", unsupported_second_sample)
    output = tmp_path / "results"
    with pytest.raises(FloatingPointError, match="model support"):
        run_analysis(
            FIXTURE / "samples.yaml",
            prepared_data,
            FIXTURE / "inputs/target",
            output,
            "paired",
        )
    assert len(processed) == 2
    assert not output.exists()
    assert not list(tmp_path.glob(".excavator2-*"))


def test_manifest_checksums_and_target_identity(tmp_path, prepared_data):
    prepared = tmp_path / "prepared"
    shutil.copytree(FIXTURE / "inputs/prepared", prepared)
    manifest = json.loads((prepared / "manifest.json").read_text())
    manifest["target_id"] = "wrong-target"
    (prepared / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="identities"):
        run_analysis(
            FIXTURE / "samples.yaml",
            prepared,
            FIXTURE / "inputs/target",
            tmp_path / "output",
            "paired",
        )
    (prepared / next(iter(manifest["files"]))).write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="checksum"):
        load_manifest(prepared, "prepared")


def test_nested_windows_are_rejected_before_analysis(tmp_path):
    prepared = tmp_path / "prepared"
    shutil.copytree(FIXTURE / "inputs/prepared", prepared)
    manifest = json.loads((prepared / "manifest.json").read_text())
    for filename in manifest["files"]:
        with np.load(prepared / filename, allow_pickle=False) as data:
            arrays = dict(data)
        # Original geometry: a wide window enclosing its successor.
        arrays["matrix"][0, 3] = arrays["matrix"][1, 3]
        np.savez_compressed(prepared / filename, **arrays)
        manifest["files"][filename] = digest(prepared / filename)
    (prepared / "manifest.json").write_text(json.dumps(manifest))
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="overlap on chr1; regenerate the target"):
        run_analysis(
            FIXTURE / "samples.yaml", prepared, FIXTURE / "inputs/target", output, "paired"
        )
    assert not output.exists()


def test_design_order_and_duplicate_yaml(tmp_path):
    samples = {"T2": "test2", "C2": "control2", "C1": "control1", "T1": "test1"}
    assert experimental_design(samples, "paired") == [
        ("test1", ["control1"]),
        ("test2", ["control2"]),
    ]
    assert experimental_design(samples, "pooled") == [
        ("test2", ["control2", "control1"]),
        ("test1", ["control2", "control1"]),
    ]
    with pytest.raises(ValueError, match="matching"):
        experimental_design({"C2": "control", "T1": "test"}, "paired")
    path = tmp_path / "samples.yaml"
    path.write_text("T1: one\nT1: two\n")
    with pytest.raises(ValueError, match="duplicate"):
        read_yaml(path)


def test_summary_preserves_cross_chromosome_boundary():
    rows = np.array(
        [["chr1", "1", "1", "2", "0", "0", "IN"], ["chr2", "3", "3", "4", "0", "0", "OUT"]]
    )
    starts, ends, values = summarize(rows, [0, 0])
    assert_array_equal(starts, [0, 1])
    assert_array_equal(ends, [1, 2])
    assert_array_equal(values, [0, 0])


def test_cli_records_per_sample_r_state_without_changing_unique_calls(tmp_path, prepared_data):
    with np.load(FIXTURE.parent / "legacy-fastcall/edges/ties.npz") as data:
        seed = data["seed_before"]
    seeds = tmp_path / "seeds.json"
    seeds.write_text(json.dumps({name: seed.tolist() for name in ("Test1", "Test2")}))
    output = tmp_path / "results"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "excavator2",
            "analyze",
            "--samples",
            str(FIXTURE / "samples.yaml"),
            "--input",
            str(prepared_data),
            "--target",
            str(FIXTURE / "inputs/target"),
            "--output",
            str(output),
            "--experiment",
            "paired",
            "--r-seed-states",
            str(seeds),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    for sample in ("Test1", "Test2"):
        with np.load(output / "Results" / sample / "checkpoints.npz") as data:
            assert_array_equal(data["r_seed_before"], seed)
            assert data["r_seed_after"].shape == (626,)
    manifest = json.loads((output / "manifest.json").read_text())
    from excavator2.artifacts import digest

    assert manifest["r_seed_states_sha256"] == digest(seeds)


def test_analysis_replays_original_ties_independently_per_sample(
    tmp_path, monkeypatch, prepared_data
):
    from dataclasses import replace

    from excavator2 import analyze

    with np.load(FIXTURE.parent / "legacy-fastcall/edges/ties.npz") as archive:
        expected = {key: archive[key] for key in archive.files}
    real_fit = analyze.fit_fastcall

    def tied_fit(*args, **kwargs):
        fit = real_fit(*args, **kwargs)
        assert len(fit.posterior) >= 4
        posterior = np.zeros_like(fit.posterior)
        posterior[:, 2] = 1
        posterior[:4] = expected["posterior"][:4]
        return replace(fit, posterior=posterior)

    monkeypatch.setattr(analyze, "fit_fastcall", tied_fit)
    seeds = tmp_path / "seeds.json"
    seeds.write_text(
        json.dumps({name: expected["seed_before"].tolist() for name in ("Test1", "Test2")})
    )
    output = tmp_path / "results"
    run_analysis(
        FIXTURE / "samples.yaml",
        prepared_data,
        FIXTURE / "inputs/target",
        output,
        "paired",
        r_seed_states=seeds,
    )
    for sample in ("Test1", "Test2"):
        with np.load(output / "Results" / sample / "checkpoints.npz") as data:
            assert_array_equal(data["labels"][:4], expected["calls"][:4, 0])
            from excavator2.fastcall import assign_labels

            replay = assign_labels(data["posterior"], r_seed=expected["seed_before"])
            assert_array_equal(data["labels"], replay.labels)
            assert_array_equal(data["r_seed_after"], replay.r_seed)


@pytest.mark.parametrize("states", [{}, {"Unknown": []}, {"Test1": [], "Test2": []}])
def test_invalid_r_states_do_not_publish_results(tmp_path, states, prepared_data):
    seeds = tmp_path / "seeds.json"
    seeds.write_text(json.dumps(states))
    output = tmp_path / "results"
    with pytest.raises(ValueError, match="R"):
        run_analysis(
            FIXTURE / "samples.yaml",
            prepared_data,
            FIXTURE / "inputs/target",
            output,
            "paired",
            r_seed_states=seeds,
        )
    assert not output.exists()


def test_legacy_converter_preserves_character_values(tmp_path):
    from excavator2.artifacts import convert_legacy

    def export(matrix, folder):
        folder.mkdir(parents=True)
        (folder / "metadata.json").write_text(
            json.dumps(
                {
                    "type": "character",
                    "length": matrix.size,
                    "dim": list(matrix.shape),
                    "order": "F",
                    "endian": "little",
                }
            )
        )
        (folder / "values.json").write_text(json.dumps(matrix.ravel(order="F").tolist()))
        np.zeros(matrix.size, dtype="<i4").tofile(folder / "missing.bin")

    with np.load(FIXTURE / "inputs/prepared/sample-0.npz") as f:
        matrix = f["matrix"]
    export(matrix, tmp_path / "exports/prepared/Control/RCNorm/Control.NRC.RData/MatrixNorm")
    for index, chromosome in enumerate(["chr1", "chr2"]):
        with np.load(FIXTURE / f"inputs/target/reference-{index}.npz") as f:
            export(f["matrix"], tmp_path / f"exports/target/FRB/FRB.{chromosome}.RData/FRBData")
    (tmp_path / "chromosomes.txt").write_text("chr1 chr2")
    (tmp_path / "centromeres.txt").write_text(
        "chrom\tstart\tend\nchr1\t7000\t11000\nchr2\t7000\t11000\n"
    )
    output = tmp_path / "converted"
    convert_legacy(
        tmp_path / "exports/prepared",
        tmp_path / "exports/target",
        tmp_path / "chromosomes.txt",
        tmp_path / "centromeres.txt",
        "synthetic",
        output,
    )
    prepared = load_manifest(output / "prepared", "prepared")
    target = load_manifest(output / "target", "target")
    assert prepared["target_id"] == target["target_id"]
    with np.load(output / "prepared/sample-0.npz") as f:
        assert_array_equal(f["matrix"], matrix)


@pytest.mark.parametrize("purity", [0.5, 0.05])
def test_analysis_passes_corrected_values_to_all_outputs(tmp_path, prepared_data, purity):
    from excavator2.analyze import DEFAULTS
    from excavator2.fastcall import correct_cellularity

    parameters = tmp_path / "parameters.json"
    params = {key: dict(value) for key, value in DEFAULTS.items()}
    params["FastCall"]["Cellularity"] = purity
    parameters.write_text(json.dumps(params))
    output = tmp_path / "corrected"
    run_analysis(
        FIXTURE / "samples.yaml",
        prepared_data,
        FIXTURE / "inputs/target",
        output,
        "paired",
        parameters,
    )
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["absolute_copy_number"] is False
    for sample in ("Test1", "Test2"):
        folder = output / "Results" / sample
        with np.load(folder / "checkpoints.npz") as data:
            corrected = correct_cellularity(data["mixture_segments"], purity)
            assert_allclose(data["corrected_segments"], corrected)
            selected = data["labels"] != 0
            table = np.genfromtxt(
                folder / f"FastCallResults_{sample}.txt", names=True, dtype=None, encoding="utf-8"
            )
            assert_allclose(table["DECNF"], 2 * np.exp2(corrected[selected]))
            assert_allclose(table["CorrectedSegment"], corrected[selected])


def test_analysis_recenters_segments_on_ploidy(tmp_path, prepared_data):
    from excavator2.analyze import DEFAULTS

    parameters = tmp_path / "parameters.json"
    params = {key: dict(value) for key, value in DEFAULTS.items()}
    params["FastCall"]["Ploidy"] = 3
    parameters.write_text(json.dumps(params))
    for name, file in [("diploid", None), ("triploid", parameters)]:
        run_analysis(
            FIXTURE / "samples.yaml",
            prepared_data,
            FIXTURE / "inputs/target",
            tmp_path / name,
            "paired",
            file,
        )
    manifest = json.loads((tmp_path / "triploid/manifest.json").read_text())
    assert manifest["parameters"]["FastCall"]["Ploidy"] == 3
    for sample in ("Test1", "Test2"):
        checkpoints = [
            np.load(tmp_path / name / "Results" / sample / "checkpoints.npz")
            for name in ("diploid", "triploid")
        ]
        diploid, triploid = checkpoints
        # Segmentation is unchanged; only segment levels move.
        assert_array_equal(diploid["segment_ids"], triploid["segment_ids"])
        assert_allclose(triploid["mixture_segments"], diploid["mixture_segments"] + np.log2(1.5))
        for data in checkpoints:
            data.close()


def test_omitted_ploidy_defaults_to_diploid(tmp_path, prepared_data):
    parameters = tmp_path / "parameters.yaml"
    parameters.write_text(
        "HSLM: {Omega: 0.1, Theta: 1e-5, D_norm: 1e6}\n"
        "FastCall: {Cellularity: 1, d: 0.5, u: 0.35, minExons: 4}\n"
    )
    output = tmp_path / "results"
    run_analysis(
        FIXTURE / "samples.yaml",
        prepared_data,
        FIXTURE / "inputs/target",
        output,
        "paired",
        parameters,
    )
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["parameters"]["FastCall"]["Ploidy"] == 2
    assert manifest["parameters"]["FastCall"]["MergeK"] == 2


def test_events_group_segments_in_region_vcf_only(tmp_path, prepared_data):
    from excavator2.analyze import DEFAULTS

    parameters = tmp_path / "parameters.json"
    params = {key: dict(value) for key, value in DEFAULTS.items()}
    params["FastCall"]["MergeK"] = 1e9
    parameters.write_text(json.dumps(params))
    for name, file in [("default", None), ("merged", parameters)]:
        run_analysis(
            FIXTURE / "samples.yaml",
            prepared_data,
            FIXTURE / "inputs/target",
            tmp_path / name,
            "paired",
            file,
        )
    for sample in ("Test1", "Test2"):
        default, merged = (tmp_path / name / "Results" / sample for name in ("default", "merged"))
        with np.load(default / "checkpoints.npz") as data:
            # Noiseless fixture steps are far beyond two standard errors.
            assert_array_equal(data["event_ids"], np.arange(len(data["labels"])))
        with np.load(merged / "checkpoints.npz") as data:
            labels, ids = data["labels"], data["event_ids"]
            # A huge threshold merges every same-call run; fixture runs lie within one arm.
            first = np.r_[True, ids[1:] != ids[:-1]]
            assert_array_equal(labels[first], data["event_labels"])
            assert (np.diff(data["event_labels"]) != 0).all()
            called = np.flatnonzero(data["event_labels"] != 0)
            table = np.genfromtxt(
                merged / f"FastCallResults_{sample}.txt", names=True, dtype=None, encoding="utf-8"
            )
            assert_array_equal(table["EventId"], ids[labels != 0])
            assert len(np.unique(table["EventId"])) < len(table)
            vcf = (merged / f"EXCAVATORRegionCall_{sample}.vcf").read_text().splitlines()
            records = [line.split("\t") for line in vcf if not line.startswith("#")]
            assert len(records) == len(called)
            probabilities = data["event_probabilities"][called]
            fcp = probabilities[np.arange(len(called)), data["event_labels"][called] + 2]
            fields = [dict(zip(r[8].split(":"), r[9].split(":"), strict=True)) for r in records]
            assert_allclose([float(f["FCP"]) for f in fields], np.round(fcp, 2))
        # Segment outputs are unchanged by merging.
        for output in ["HSLMResults", "EXCAVATORWindowCall"]:
            suffix = "txt" if output == "HSLMResults" else "vcf"
            old = (default / f"{output}_{sample}.{suffix}").read_text().splitlines()
            new = (merged / f"{output}_{sample}.{suffix}").read_text().splitlines()
            assert [x for x in old if "fileDate" not in x] == [
                x for x in new if "fileDate" not in x
            ]


def test_analysis_help_has_no_calibration_option():
    result = subprocess.run(
        [sys.executable, "-m", "excavator2", "analyze", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "--calibration" not in result.stdout
