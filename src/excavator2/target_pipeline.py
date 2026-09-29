"""Target command orchestration and version-one artifacts; policy stays Python."""

import json
import shutil
import tempfile
from pathlib import Path

import numpy as np

from . import __version__
from .artifacts import digest, metadata_identity, read_yaml, window_metadata
from .target import target_geometry
from .target_features import target_features


def run_target(config, output, force=False):
    output = Path(output)
    if force or output.exists():
        raise ValueError(
            "target generation requires a new output directory; --force is not supported"
        )
    settings = read_yaml(config)
    reference, target_settings = settings.get("Reference"), settings.get("Target")
    if not isinstance(reference, dict) or not isinstance(target_settings, dict):
        raise ValueError("config requires Reference and Target mappings")
    assembly, name = reference.get("Assembly"), target_settings.get("Name")
    if not isinstance(assembly, str) or not assembly.strip():
        raise ValueError("Reference.Assembly must be a nonempty string")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("Target.Name must be a nonempty string")
    paths = {}
    for key in ["FASTA", "BigWig", "Chromosomes", "Centromeres", "Gaps"]:
        value = reference.get(key)
        if not isinstance(value, str):
            raise ValueError(f"Reference.{key} must be a path")
        paths[key] = Path(value).resolve(strict=True)
    if not isinstance(target_settings.get("BED"), str):
        raise ValueError("Target.BED must be a path")
    paths["BED"] = Path(target_settings["BED"]).resolve(strict=True)
    # As with the original CLI, relative reference paths use the working directory.
    target = target_geometry(
        paths["BED"], paths["Chromosomes"], paths["Gaps"], target_settings.get("Window")
    )
    chromosomes = list(dict.fromkeys(target[:, 0]))
    centro = np.loadtxt(paths["Centromeres"], dtype=str, skiprows=1, ndmin=2)
    centers = {
        row[0] if len(chromosomes[0]) >= 4 else row[0][3:]: [int(row[1]), int(row[2])]
        for row in centro
    }
    if any(chrom not in centers for chrom in chromosomes):
        raise ValueError("centromere table is missing target chromosomes")
    features = target_features(target, paths["FASTA"], paths["BigWig"])
    # A malformed legacy target may contain skipped or duplicated feature rows.
    # Preserve extraction behavior, but do not publish a falsely usable artifact.
    for chrom in chromosomes:
        windows = target[target[:, 0] == chrom]
        values = features[chrom]
        if any(len(values[key]) != len(windows) for key in ["gc", "mappability", "first_base"]):
            raise ValueError(f"legacy feature rows do not align with target windows for {chrom}")
        if not np.array_equal(values["first_base"][:, 0], windows[:, 1]):
            raise ValueError(f"legacy reference coordinates do not align for {chrom}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".excavator2-target-", dir=output.parent))
    try:
        manifest = {
            "schema": 1,
            "kind": "target",
            "target_id": metadata_identity(window_metadata(target)),
            "assembly": assembly,
            "name": name.strip(),
            "chromosomes": chromosomes,
            "centromeres": centers,
            "references": {},
            "files": {},
            "preparation": "preparation.npz",
            "window": target_settings["Window"],
            "package_version": __version__,
            "compatibility": "legacy",
            "config_sha256": digest(config),
            "inputs": {
                key: {"path": str(path), "sha256": digest(path)} for key, path in paths.items()
            },
        }
        np.savez_compressed(
            temporary / "preparation.npz",
            target=target,
            gc=np.concatenate([features[c]["gc"] for c in chromosomes]),
            mappability=np.concatenate([features[c]["mappability"] for c in chromosomes]),
        )
        manifest["files"]["preparation.npz"] = digest(temporary / "preparation.npz")
        for index, chrom in enumerate(chromosomes):
            filename = f"reference-{index}.npz"
            np.savez_compressed(temporary / filename, matrix=features[chrom]["first_base"])
            manifest["references"][chrom] = filename
            manifest["files"][filename] = digest(temporary / filename)
        np.savetxt(temporary / "Filtered.txt", target, fmt="%s", delimiter="\t")
        manifest["files"]["Filtered.txt"] = digest(temporary / "Filtered.txt")
        (temporary / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        temporary.rename(output)
    except BaseException:
        shutil.rmtree(temporary)
        raise
