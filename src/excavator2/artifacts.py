"""Minimal versioned analysis artifacts converted from the binary legacy export.

R is used once by tools/oracle/export.R. Reading the resulting artifacts and
running analysis require only Python and the compiled kernels.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
import yaml


class UniqueLoader(yaml.SafeLoader):
    pass


def _mapping(loader, node):
    pairs = loader.construct_pairs(node, deep=True)
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate YAML key: {key}")
        result[key] = value
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def read_yaml(path):
    with Path(path).open() as handle:
        value = yaml.load(handle, Loader=UniqueLoader)
    if not isinstance(value, dict):
        raise ValueError(f"expected a YAML mapping: {path}")
    return value


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_export(folder):
    """Read one atomic object exported by export.R; preserve character precision."""
    folder = Path(folder)
    meta = json.loads((folder / "metadata.json").read_text())
    if meta["order"] != "F" or meta["endian"] != "little":
        raise ValueError("unsupported legacy export layout")
    if np.fromfile(folder / "missing.bin", dtype="<i4").any():
        raise ValueError(f"missing values in legacy artifact: {folder}")
    if meta["type"] == "character":
        values = np.array(json.loads((folder / "values.json").read_text()), dtype=str)
    elif meta["type"] in ("double", "integer", "logical"):
        values = np.fromfile(
            folder / "data.bin", dtype="<f8" if meta["type"] == "double" else "<i4"
        )
    else:
        raise ValueError("unsupported legacy object type")
    if values.size != meta["length"]:
        raise ValueError("truncated legacy export")
    return values.reshape(meta["dim"], order="F") if meta["dim"] is not None else values


def load_manifest(folder, kind):
    folder = Path(folder)
    manifest = json.loads((folder / "manifest.json").read_text())
    if manifest.get("schema") != 1 or manifest.get("kind") != kind:
        raise ValueError(f"expected version 1 {kind} artifact")
    for filename, checksum in manifest["files"].items():
        if Path(filename).name != filename or digest(folder / filename) != checksum:
            raise ValueError(f"invalid artifact file or checksum: {filename}")
    return manifest


def convert_legacy(prepared_exports, target_exports, chromosomes, centromeres, assembly, output):
    """Convert already exported NRC/FRB objects; never reinterpret raw RData."""
    output = Path(output)
    if output.exists():
        raise ValueError("conversion output must not exist")
    matrices = {}
    identity = None
    for path in sorted(Path(prepared_exports).glob("*/RCNorm/*.NRC.RData/MatrixNorm")):
        sample = path.parents[2].name
        matrix = read_export(path)
        if matrix.ndim != 2 or matrix.shape[1] != 7:
            raise ValueError("expected a seven-column MatrixNorm")
        # Canonical per-column strings avoid a hash depending on NumPy string width.
        current = hashlib.sha256(
            json.dumps(matrix[:, [0, 1, 2, 3, 4, 6]].tolist(), separators=(",", ":")).encode()
        ).hexdigest()
        if identity is not None and current != identity:
            raise ValueError("prepared samples have mismatched target windows")
        identity = current
        matrices[sample] = matrix
    if not matrices:
        raise ValueError("no exported MatrixNorm objects found")
    chroms = Path(chromosomes).read_text().split()
    if len(set(chroms)) != len(chroms):
        raise ValueError("duplicate target chromosomes")
    centro = np.loadtxt(centromeres, dtype=str, skiprows=1, ndmin=2)
    centers = {
        row[0] if len(chroms[0]) >= 4 else row[0][3:]: [int(row[1]), int(row[2])] for row in centro
    }
    refs = {}
    for chromosome in chroms:
        if chromosome not in centers:
            raise ValueError(f"missing centromere: {chromosome}")
        refs[chromosome] = read_export(
            Path(target_exports) / "FRB" / f"FRB.{chromosome}.RData" / "FRBData"
        )
    output.mkdir(parents=True)
    for kind in ["prepared", "target"]:
        (output / kind).mkdir()
    prepared = {"schema": 1, "kind": "prepared", "target_id": identity, "samples": {}, "files": {}}
    for index, (sample, matrix) in enumerate(matrices.items()):
        filename = f"sample-{index}.npz"
        np.savez_compressed(output / "prepared" / filename, matrix=matrix)
        prepared["samples"][sample] = filename
        prepared["files"][filename] = digest(output / "prepared" / filename)
    target = {
        "schema": 1,
        "kind": "target",
        "target_id": identity,
        "assembly": assembly,
        "chromosomes": chroms,
        "centromeres": centers,
        "references": {},
        "files": {},
    }
    for index, (chromosome, matrix) in enumerate(refs.items()):
        filename = f"reference-{index}.npz"
        np.savez_compressed(output / "target" / filename, matrix=matrix)
        target["references"][chromosome] = filename
        target["files"][filename] = digest(output / "target" / filename)
    for kind, manifest in [("prepared", prepared), ("target", target)]:
        (output / kind / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
