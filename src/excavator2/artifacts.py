"""Minimal versioned analysis artifacts converted from the binary legacy export.

R is used once by tools/oracle/export.R. Reading the resulting artifacts and
running analysis require only Python and the compiled kernels.
"""

import hashlib
import json
from decimal import Decimal
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
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


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
    exported_samples = (
        sorted(Path(prepared_exports).glob("*/RCNorm/*.NRC.RData/MatrixNorm"))
        if prepared_exports is not None
        else []
    )
    for path in exported_samples:
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
    features = load_preparation_exports(target_exports, chroms)
    if features is not None:
        metadata = window_metadata(features["target"])
        current = metadata_identity(metadata)
        if identity is not None and current != identity:
            raise ValueError("target features do not match prepared window metadata")
        identity = current
    if identity is None:
        raise ValueError("no prepared samples or complete target feature exports found")
    output.mkdir(parents=True)
    for kind in ["prepared", "target"] if matrices else ["target"]:
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
    if features is not None:
        filename = "preparation.npz"
        np.savez_compressed(output / "target" / filename, **features)
        target["preparation"] = filename
        target["files"][filename] = digest(output / "target" / filename)
    manifests = [("prepared", prepared), ("target", target)] if matrices else [("target", target)]
    for kind, manifest in manifests:
        (output / kind / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def fixed_number(value):
    """R's 15 significant digits with scipen=20 at the normalization boundary."""
    return format(Decimal(format(float(value), ".15g")), "f")


def window_metadata(target):
    start, end = coordinate_values(target[:, 1]), coordinate_values(target[:, 2])
    return np.column_stack(
        [
            target[:, 0],
            [fixed_number(x) for x in (start + end) / 2],
            start.astype(str),
            end.astype(str),
            target[:, 3],
            target[:, 4],
        ]
    )


def metadata_identity(metadata):
    return hashlib.sha256(json.dumps(metadata.tolist(), separators=(",", ":")).encode()).hexdigest()


def load_preparation_exports(folder, chromosomes):
    folder = Path(folder)
    targets = list(folder.glob("*.RData/MyTarget"))
    if not targets:
        return None  # Existing M3-only converted inputs remain readable.
    if len(targets) != 1:
        raise ValueError("expected one exported MyTarget")
    target = read_export(targets[0])
    gcfiles = sorted((folder / "GCC").glob("*.RData"))
    mapfiles = sorted((folder / "MAP").glob("*.RData"))
    gc, maps = [], []
    for chromosome in chromosomes:
        gi = [f.name for f in gcfiles].index(f"GCC.{chromosome}.RData")
        mi = [f.name for f in mapfiles].index(f"Map.{chromosome}.RData")
        # Preserve loadTarget's crossed indices, including its directory ordering.
        gc.append(read_export(gcfiles[mi] / "GCContent"))
        maps.append(read_export(mapfiles[gi] / "MapMed"))
    return {"target": target, "gc": np.concatenate(gc), "mappability": np.concatenate(maps)}


def coordinate_values(values):
    """R as.integer character parsing, with defined signed-32-bit input bounds."""
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all() or (np.abs(values) > np.iinfo(np.int32).max).any():
        raise ValueError("coordinates exceed the finite legacy integer range")
    return np.trunc(values).astype(np.int64)
