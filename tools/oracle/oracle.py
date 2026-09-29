"""Run and compare the pinned legacy fixture. Never changes legacy scientific code."""

import argparse
import hashlib
import json
import shutil
import subprocess
import tarfile
import tempfile
from datetime import UTC, datetime
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
BASELINE = json.loads((Path(__file__).with_name("baseline.json")).read_text())


def file_record(path: Path) -> dict:
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"size": path.stat().st_size, "sha256": digest}


def save_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def container_command(output: Path) -> list[str]:
    return [
        "docker",
        "run",
        "--rm",
        "--platform",
        BASELINE["platform"],
        "--network",
        "none",
        "--entrypoint",
        "micromamba",
        "--mount",
        f"type=bind,src={ROOT},dst=/repo,readonly",
        "--mount",
        f"type=bind,src={output},dst=/work",
        BASELINE["image"],
        "run",
        "-n",
        "excavator2",
    ]


def run(output: Path, threads: int) -> int:
    config_path = ROOT / BASELINE["config"]
    samples_path = ROOT / BASELINE["samples"]
    design_path = ROOT / BASELINE["design"]
    config = yaml.safe_load(config_path.read_text())
    samples = yaml.safe_load(samples_path.read_text())
    design = yaml.safe_load(design_path.read_text())
    inputs = [config_path, samples_path, design_path, ROOT / "tools/oracle/baseline.json"]
    for key in ("FASTA", "BigWig", "Chromosomes", "Centromeres", "Gaps"):
        path = ROOT / config["Reference"][key]
        inputs.append(path)
        config["Reference"][key] = "/repo/" + path.relative_to(ROOT).as_posix()
    bed = ROOT / config["Target"]["BED"]
    inputs.append(bed)
    config["Target"]["BED"] = "/repo/" + bed.relative_to(ROOT).as_posix()
    for label, filename in samples.items():
        path = ROOT / filename
        inputs.extend([path, Path(str(path) + ".bai")])
        samples[label] = "/repo/" + path.relative_to(ROOT).as_posix()
    fasta = ROOT / config["Reference"]["FASTA"].removeprefix("/repo/")
    inputs.append(Path(str(fasta) + ".fai"))
    missing = [str(path.relative_to(ROOT)) for path in inputs if not path.is_file()]
    if missing:
        raise ValueError("Missing required inputs: " + ", ".join(missing))
    subprocess.run(
        ["docker", "image", "inspect", BASELINE["image"]], check=True, stdout=subprocess.DEVNULL
    )
    records = {str(p.relative_to(ROOT)): file_record(p) for p in inputs}
    expected = json.loads(Path(__file__).with_name("inputs.sha256.json").read_text())
    if records != expected:
        raise ValueError(
            "Inputs differ from the pinned fixture inventory; do not overwrite goldens"
        )
    output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "baseline": BASELINE,
        "threads": threads,
        "status": "preparing",
        "started_utc": datetime.now(UTC).isoformat(),
        "inputs": records,
        "harness": {
            name: file_record(Path(__file__).with_name(name))
            for name in ("oracle.py", "run.sh", "export.R")
        },
        "seed_policy": "unmodified legacy RNG; no injected seed",
    }
    manifest_path = output / "manifest.json"
    save_json(manifest_path, manifest)
    try:
        (output / "harness").mkdir()
        for name in ("oracle.py", "run.sh", "export.R", "baseline.json"):
            shutil.copy2(Path(__file__).with_name(name), output / "harness" / name)
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
        manifest["source_files"] = {
            str(p.relative_to(source)): file_record(p)
            for p in sorted(source.rglob("*"))
            if p.is_file()
        }
        for filename, data in (
            ("config.yaml", config),
            ("samples.yaml", samples),
            ("design.yaml", design),
        ):
            (output / filename).write_text(yaml.safe_dump(data, sort_keys=False))
        target = (
            f"/work/target/{config['Reference']['Assembly']}/"
            f"{config['Target']['Name']}/w_{config['Target']['Window']}"
        )
        (output / "target-path.txt").write_text(target + "\n")
        command = container_command(output) + ["bash", "/work/harness/run.sh", str(threads)]
        manifest.update(status="running", command=command)
        save_json(manifest_path, manifest)
        with (output / "container.log").open("w") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        manifest["returncode"] = result.returncode
        if result.returncode:
            raise RuntimeError("Legacy execution failed; inspect container/stage/compile logs")
        # Legacy parent scripts do not reliably propagate child failures.
        required = [output / target.removeprefix("/work/") / "Filtered.txt"]
        required += [
            output / "prepared" / label / "RCNorm" / f"{label}.NRC.RData" for label in samples
        ]
        for tag, label in design.items():
            if tag.startswith("T"):
                folder = output / "results/Results" / label
                required += [
                    folder / f"{prefix}_{label}.{suffix}"
                    for prefix, suffix in (
                        ("HSLMResults", "txt"),
                        ("FastCallResults", "txt"),
                        ("EXCAVATORWindowCall", "vcf"),
                        ("EXCAVATORRegionCall", "vcf"),
                    )
                ]
        absent = [
            str(p.relative_to(output)) for p in required if not p.is_file() or p.stat().st_size == 0
        ]
        if absent:
            raise RuntimeError("Missing/empty legacy outputs: " + ", ".join(absent))
        with (output / "export.log").open("w") as log:
            subprocess.run(
                container_command(output)
                + ["Rscript", "/work/harness/export.R", "/work", "/work/checkpoints"],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        manifest["artifacts"] = {
            str(p.relative_to(output)): file_record(p)
            for folder in ("target", "prepared", "results", "checkpoints", "runtime")
            for p in sorted((output / folder).rglob("*"))
            if p.is_file()
        }
        manifest["status"] = "complete"
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        manifest.update(status="failed", error=str(error))
    finally:
        manifest["finished_utc"] = datetime.now(UTC).isoformat()
        save_json(manifest_path, manifest)
    print(json.dumps({"status": manifest["status"], "manifest": str(manifest_path)}))
    return 0 if manifest["status"] == "complete" else 1


def comparable_files(root: Path) -> dict[str, bytes]:
    """Compare every exported object and final call table in original row order."""
    paths = list((root / "checkpoints").rglob("*"))
    paths += list((root / "results/Results").rglob("*"))
    result = {}
    for path in sorted(paths):
        if not path.is_file() or path.suffix not in (".bin", ".json", ".txt", ".vcf"):
            continue
        data = path.read_bytes()
        if path.suffix == ".vcf":
            data = b"\n".join(
                line for line in data.split(b"\n") if not line.startswith(b"##fileDate=")
            )
        result[str(path.relative_to(root))] = data
    return result


def compare(left: Path, right: Path, report: Path) -> int:
    manifests = [json.loads((root / "manifest.json").read_text()) for root in (left, right)]
    for root, manifest in zip((left, right), manifests, strict=True):
        if manifest["status"] != "complete":
            raise ValueError(f"Cannot compare incomplete run: {root}")
    for field in ("baseline", "inputs", "source_files"):
        if field not in manifests[0] or field not in manifests[1]:
            raise ValueError(f"Missing comparison provenance: {field}")
        if manifests[0][field] != manifests[1][field]:
            raise ValueError(f"Different comparison provenance: {field}")
    for name in ("run.sh", "export.R"):
        if manifests[0].get("harness", {}).get(name) != manifests[1].get("harness", {}).get(name):
            raise ValueError(f"Different execution/export harness: {name}")
    a, b = comparable_files(left), comparable_files(right)
    if not a or not b:
        raise ValueError("No comparison artifacts")
    differences = [key for key in sorted(a.keys() | b.keys()) if a.get(key) != b.get(key)]
    result = {
        "left": str(left),
        "right": str(right),
        "compared_files": len(a.keys() | b.keys()),
        "differences": differences,
        "equal": not differences,
        "ignored_fields": ["VCF ##fileDate"],
        "scope": "exported RData objects and final txt/vcf; plots/logs excluded",
    }
    save_json(report, result)
    print(json.dumps(result, indent=2))
    return int(bool(differences))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run_parser = commands.add_parser("run")
    run_parser.add_argument("--output", required=True, type=Path)
    run_parser.add_argument("--threads", type=int, default=1)
    compare_parser = commands.add_parser("compare")
    compare_parser.add_argument("left", type=Path)
    compare_parser.add_argument("right", type=Path)
    compare_parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    try:
        if args.command == "run":
            if args.threads < 1:
                parser.error("--threads must be positive")
            return run(args.output.resolve(), args.threads)
        return compare(args.left.resolve(), args.right.resolve(), args.report)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"oracle: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
