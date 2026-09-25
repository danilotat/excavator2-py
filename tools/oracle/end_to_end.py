"""Fresh-reference, non-normal target/prepare/analyze concordance for both designs."""

import argparse
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pyBigWig
import pysam
import yaml
from compare_analysis import compare as compare_analysis
from compare_target import compare as compare_target
from oracle import BASELINE, ROOT, container_command, file_record, save_json

from excavator2.artifacts import convert_legacy, load_manifest, read_export


def make_inputs(output):
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    inputs = output / "inputs"
    inputs.mkdir()
    (inputs / "reference.fa").write_text(
        "".join(f">{c}\n" + "ACGT" * 6000 + "\n" for c in chromosomes)
    )
    pysam.faidx(str(inputs / "reference.fa"))
    with pyBigWig.open(str(inputs / "reference.bw"), "w") as bw:
        bw.addHeader([(c, 24000) for c in chromosomes])
        bw.addEntries(chromosomes, [0] * 23, ends=[24000] * 23, values=[1.0] * 23)
    (inputs / "coordinates.tsv").write_text("".join(f"{c}\t0\t20000\n" for c in chromosomes))
    (inputs / "centromeres.tsv").write_text(
        "chrom\tstart\tend\n" + "".join(f"{c}\t9000\t11000\n" for c in chromosomes)
    )
    (inputs / "gaps.tsv").write_text(
        "id\tchrom\tstart\tend\n"
        + "".join(f"0\t{c}\t0\t1\n0\t{c}\t9000\t11000\n" for c in chromosomes)
        + "0\tchr1_alt\t0\t1\n"
    )
    (inputs / "target.bed").write_text(
        "".join(f"chr1\t{s}\t{s + 40}\n" for s in range(500, 19501, 1000))
    )
    references = {
        "FASTA": "reference.fa",
        "BigWig": "reference.bw",
        "Chromosomes": "coordinates.tsv",
        "Centromeres": "centromeres.tsv",
        "Gaps": "gaps.tsv",
    }
    for kind, base in [("current", str(inputs)), ("legacy", "/work/inputs")]:
        config = {
            "Reference": {
                "Assembly": "synthetic",
                **{k: f"{base}/{v}" for k, v in references.items()},
            },
            "Target": {"Name": "panel", "BED": f"{base}/target.bed", "Window": 100},
        }
        (output / f"{kind}-config.yaml").write_text(yaml.safe_dump(config))
    samples = {}
    for sample in ["Control1", "Control2", "Test1"]:
        path = inputs / f"{sample}.bam"
        header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": c, "LN": 24000} for c in chromosomes]}
        with pysam.AlignmentFile(path, "wb", header=header) as bam:
            for index, chrom in enumerate(chromosomes):
                for position in range(25, 22000, 100):
                    count = 30 if sample == "Control2" else 25
                    if sample == "Test1" and chrom in ["chr1", "chr2"]:
                        if 1000 <= position < 5000:
                            count = 10
                        elif 12000 <= position < 16000:
                            count = 50
                    for j in range(count):
                        read = pysam.AlignedSegment()
                        read.query_name = f"{index}_{position}_{j}"
                        read.reference_id = index
                        read.reference_start = position
                        read.mapping_quality = 60
                        read.query_sequence = "A"
                        read.cigarstring = "1M"
                        bam.write(read)
        pysam.index(str(path))
        samples[sample] = str(path)
    (output / "current-samples.yaml").write_text(yaml.safe_dump(samples))
    (output / "legacy-samples.yaml").write_text(
        yaml.safe_dump({k: f"/work/inputs/{k}.bam" for k in samples})
    )
    (output / "paired.yaml").write_text("T1: Test1\nC1: Control1\n")
    (output / "pooling.yaml").write_text("T1: Test1\nC1: Control1\nC2: Control2\n")


def run(output):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = {
        "passed": False,
        "baseline": BASELINE,
        "scope": "scientific data and calls; plots excluded",
    }
    try:
        make_inputs(output)
        source = output / "source"
        source.mkdir()
        with tempfile.TemporaryFile() as archive:
            subprocess.run(
                ["git", "archive", BASELINE["source_commit"], "excavator2"],
                cwd=ROOT,
                stdout=archive,
                check=True,
            )
            archive.seek(0)
            with tarfile.open(fileobj=archive) as tar:
                tar.extractall(source, filter="data")
        report["inputs"] = {
            str(p.relative_to(output)): file_record(p) for p in (output / "inputs").iterdir()
        }
        report["harness"] = {
            name: file_record(Path(__file__).with_name(name))
            for name in [
                "end_to_end.py",
                "run_end_to_end.sh",
                "export.R",
                "compare_target.py",
                "compare_analysis.py",
            ]
        }
        with (output / "legacy.log").open("w") as log:
            subprocess.run(
                container_command(output) + ["bash", "/repo/tools/oracle/run_end_to_end.sh"],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        legacy_target = output / "legacy-target/synthetic/panel/w_100"
        convert_legacy(
            output / "exports/prepared",
            output / "exports/target",
            legacy_target / "panel_chromosome.txt",
            output / "inputs/centromeres.tsv",
            "synthetic",
            output / "converted",
        )
        commands = [
            [
                "target",
                "--config",
                str(output / "current-config.yaml"),
                "--output",
                str(output / "current-target"),
            ],
            [
                "prepare",
                "--samples",
                str(output / "current-samples.yaml"),
                "--target",
                str(output / "current-target"),
                "--output",
                str(output / "current-prepared"),
                "--threads",
                "2",
            ],
        ]
        commands += [
            [
                "analyze",
                "--samples",
                str(output / f"{mode}.yaml"),
                "--target",
                str(output / "current-target"),
                "--input",
                str(output / "current-prepared"),
                "--output",
                str(output / f"current-{mode}"),
                "--experiment",
                mode,
            ]
            for mode in ["paired", "pooling"]
        ]
        for index, command in enumerate(commands):
            with (output / f"current-{index}.log").open("w") as log:
                subprocess.run(
                    [sys.executable, "-I", "-m", "excavator2", *command],
                    cwd=output,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
        report["target"] = compare_target(output / "converted/target", output / "current-target")
        prepared = load_manifest(output / "current-prepared", "prepared")
        report["prepared"] = {}
        for sample, filename in prepared["samples"].items():
            original = read_export(
                output / f"exports/prepared/{sample}/RCNorm/{sample}.NRC.RData/MatrixNorm"
            )
            with np.load(output / "current-prepared" / filename) as current:
                np.testing.assert_array_equal(current["matrix"], original, err_msg=sample)
                counts = []
                for chrom in load_manifest(output / "current-target", "target")["chromosomes"]:
                    counts.append(
                        read_export(
                            output / f"exports/prepared/{sample}/RC/{sample}.RC.{chrom}.RData/RC"
                        )
                    )
                np.testing.assert_array_equal(
                    current["counts"], np.concatenate(counts), err_msg=sample
                )
                report["prepared"][sample] = {"windows": len(original), "exact": True}
        report["designs"] = {}
        for mode in ["paired", "pooling"]:
            result = compare_analysis(output / f"legacy-{mode}", output / f"current-{mode}")
            manifest = json.loads((output / f"current-{mode}/manifest.json").read_text())
            calls = manifest["samples"]["Test1"]["calls"]
            with np.load(output / f"current-{mode}/Results/Test1/checkpoints.npz") as checkpoint:
                labels = checkpoint["labels"]
            if (
                calls != 4
                or not (np.any(labels < 0) and np.any(labels > 0))
                or len(result["files"]) != 4
            ):
                raise ValueError(
                    f"{mode}: fixture must produce four loss/gain calls and all four outputs"
                )
            result["non_normal_calls"] = calls
            report["designs"][mode] = result
        report["legacy_plotting"] = {
            "qualified": False,
            "failed_modes": [
                mode
                for mode in ["paired", "pooling"]
                if "Couldn't generate plots" in (output / f"{mode}.log").read_text()
            ],
        }
        report["passed"] = True
    except Exception as error:
        report["error"] = str(error)
    finally:
        save_json(output / "concordance.json", report)
    print(json.dumps({"passed": report["passed"], "report": str(output / "concordance.json")}))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    raise SystemExit(run(parser.parse_args().output))
