"""Prepare BAM counts with original size, MAP and GC normalization."""

import json
import re
import shutil
import tempfile
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
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
from .reads import RegionCounter, region_jobs


def run_preparation(samples, target_folder, output, threads=1, mapq=20, force=False, progress=None):
    target_folder, output = Path(target_folder), Path(output)
    if force or output.exists():
        raise ValueError("preparation requires a new output directory; --force is not supported")
    if threads < 1 or not 0 <= mapq <= 255:
        raise ValueError("invalid threads or MAPQ")
    if progress:
        progress("Loading target and sample configuration")
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
    jobs = list(region_jobs(target, manifest["chromosomes"], threads))
    workers = min(threads, len(paths) * len(jobs))
    if progress:
        noun = "worker" if workers == 1 else "workers"
        progress(f"Preparing {len(paths)} samples with {workers} {noun} across target regions")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".excavator2-prepare-", dir=output.parent))
    try:

        def prepare_sample(index, name, counts):
            trace = normalize(counts, target, gc, mappability)
            normalized = [fixed_number(x) for x in trace["normalized"]]
            matrix = np.column_stack([metadata[:, :5], normalized, metadata[:, 5]])
            filename = f"sample-{index}.npz"
            np.savez_compressed(temporary / filename, matrix=matrix, counts=counts, **trace)
            return name, filename, digest(temporary / filename)

        items = list(paths.items())
        counts = {}
        remaining = [len(jobs)] * len(items)
        completed = {}
        counter = RegionCounter()
        with ThreadPoolExecutor(max_workers=workers) as pool:
            work = iter((index, path, job) for index, (_, path) in enumerate(items) for job in jobs)
            pending = {}

            def queue_jobs():
                # Bound queued work and let completed samples normalize promptly.
                while len(pending) < 2 * workers:
                    item = next(work, None)
                    if item is None:
                        break
                    index, path, job = item
                    pending[pool.submit(counter, path, job, mapq)] = (index, job[0])

            queue_jobs()
            while pending:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    index, indices = pending.pop(future)
                    name = items[index][0]
                    try:
                        result = future.result()
                    except Exception as error:
                        for queued in pending:
                            queued.cancel()
                        raise ValueError(
                            f"preparation failed for sample {name}: {error}"
                        ) from error
                    if indices is not None:
                        if index not in counts:
                            counts[index] = np.zeros(len(target), dtype=np.int64)
                        counts[index][indices] = result
                        remaining[index] -= 1
                        if remaining[index] == 0:
                            future = pool.submit(prepare_sample, index, name, counts.pop(index))
                            pending[future] = (index, None)
                    else:
                        completed[index] = result
                        if progress:
                            width = 20
                            filled = width * len(completed) // len(paths)
                            bar = "=" * filled + "-" * (width - filled)
                            progress(f"[{bar}] {len(completed)}/{len(paths)} Prepared {name}")
                queue_jobs()
        results = [completed[index] for index in range(len(items))]
        if progress:
            progress("Writing preparation manifest")
        result = {
            "schema": 1,
            "kind": "prepared",
            "depth_policy": "size-map-gc-normalized-v1",
            "target_id": manifest["target_id"],
            "samples": {name: filename for name, filename, _ in results},
            "files": {filename: checksum for _, filename, checksum in results},
            "threads": threads,
            "requested_mapq": mapq,
            "filter_policy": "flags & 1028 == 0; MAPQ >= requested_mapq",
            "count_policy": "independent-inclusive-v1",
            "count_backend": "cpp.htslib-prefix-counts",
        }
        (temporary / "manifest.json").write_text(json.dumps(result, indent=2) + "\n")
        temporary.rename(output)
        if progress:
            progress(f"Complete: {output}")
    except BaseException:
        shutil.rmtree(temporary)
        raise
