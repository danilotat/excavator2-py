"""Prepare legacy-compatible normalized BAM counts using Python and NumPy."""

import json
import re
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from .artifacts import (
    digest,
    fixed_number,
    load_manifest,
    metadata_identity,
    read_yaml,
    window_metadata,
)
from .normalization import normalize
from .reads import count_bam


def run_preparation(samples, target_folder, output, threads=1, mapq=20, force=False):
    target_folder, output = Path(target_folder), Path(output)
    if force or output.exists():
        raise ValueError("preparation requires a new output directory; --force is not supported")
    if threads < 1 or not 0 <= mapq <= 255:
        raise ValueError("invalid threads or MAPQ")
    manifest = load_manifest(target_folder, "target")
    feature_file = manifest.get("preparation")
    if feature_file not in manifest["files"]:
        raise ValueError(
            "target lacks preparation features; convert original target GCC/MAP exports"
        )
    with np.load(target_folder / feature_file, allow_pickle=False) as archive:
        target, gc, mappability = (archive[k] for k in ["target", "gc", "mappability"])
    metadata = window_metadata(target)
    if metadata_identity(metadata) != manifest["target_id"]:
        raise ValueError("preparation features have a different target identity")
    design = read_yaml(samples)
    if not design or not all(
        isinstance(k, str)
        and re.fullmatch(r"[\w.-]+", k)
        and k not in (".", "..")
        and isinstance(v, str)
        for k, v in design.items()
    ):
        raise ValueError("samples must map plain sample names to BAM paths")
    # Like the original CLI, relative paths are resolved from the working directory.
    paths = {name: Path(value).resolve(strict=True) for name, value in design.items()}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".excavator2-prepare-", dir=output.parent))
    try:

        def prepare_sample(item):
            index, (name, path) = item
            counts = count_bam(path, target, manifest["chromosomes"])
            trace = normalize(counts, target, gc, mappability)
            normalized = [fixed_number(x) for x in trace["normalized"]]
            matrix = np.column_stack([metadata[:, :5], normalized, metadata[:, 5]])
            filename = f"sample-{index}.npz"
            np.savez_compressed(temporary / filename, matrix=matrix, counts=counts, **trace)
            return name, filename, digest(temporary / filename), digest(path)

        with ThreadPoolExecutor(max_workers=min(threads, len(paths))) as pool:
            results = list(pool.map(prepare_sample, enumerate(paths.items())))
        result = {
            "schema": 1,
            "kind": "prepared",
            "target_id": manifest["target_id"],
            "samples": {name: filename for name, filename, _, _ in results},
            "files": {filename: checksum for _, filename, checksum, _ in results},
            "bam_sha256": {name: checksum for name, _, _, checksum in results},
            "target_manifest_sha256": digest(target_folder / "manifest.json"),
            "samples_sha256": digest(samples),
            "threads": threads,
            "requested_mapq": mapq,
            "filter_policy": "legacy flags & 1028 == 0; all MAPQ retained",
            "count_backend": "numpy.searchsorted",
            "chunk_size": 500000,
        }
        (temporary / "manifest.json").write_text(json.dumps(result, indent=2) + "\n")
        temporary.rename(output)
    except BaseException:
        shutil.rmtree(temporary)
        raise
