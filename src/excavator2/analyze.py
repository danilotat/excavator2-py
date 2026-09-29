"""Analysis from exported legacy prepared counts; scientific policy stays in Python."""

import json
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np

from .artifacts import digest, load_manifest, read_yaml
from .fastcall import assign_labels, correct_cellularity, fit_fastcall
from .hslm import estimate_parameters, segment
from .writers import number, write_results

DEFAULTS = {
    "HSLM": {"Omega": 0.1, "Theta": 1e-5, "D_norm": 1e6},
    "FastCall": {"Cellularity": 1.0, "d": 0.5, "u": 0.35, "minExons": 4},
}


def experimental_design(samples, experiment):
    """Paired numeric order or pooled controls in YAML insertion order."""
    if not samples or not all(
        isinstance(k, str)
        and re.fullmatch(r"[CT][0-9]+", k)
        and isinstance(v, str)
        and re.fullmatch(r"[\w.-]+", v)
        and v not in (".", "..")
        for k, v in samples.items()
    ):
        raise ValueError("samples require C<number>/T<number> labels and plain sample names")
    if len(set(samples.values())) != len(samples):
        raise ValueError("duplicate sample names")
    controls = [(int(k[1:]), v) for k, v in samples.items() if k.startswith("C")]
    tests = [(int(k[1:]), v) for k, v in samples.items() if k.startswith("T")]
    if not controls or not tests:
        raise ValueError("both controls and tests are required")
    if experiment == "paired":
        controls, tests = sorted(controls), sorted(tests)
        if [i for i, _ in controls] != [i for i, _ in tests] or len({i for i, _ in tests}) != len(
            tests
        ):
            raise ValueError("paired labels must have matching unique numeric suffixes")
        return [(test, [control]) for (_, test), (_, control) in zip(tests, controls, strict=True)]
    if experiment in ("pooling", "pooled"):
        return [(test, [control for _, control in controls]) for _, test in tests]
    raise ValueError("unsupported experimental design")


def ratios(test, controls):
    """Legacy ordered pooling, character roundtrip, then separate IN/OUT centering."""
    counts = controls[0][:, 5].astype(float)
    for control in controls[1:]:
        counts = counts + control[:, 5].astype(float)
    if len(controls) > 1:
        counts /= len(controls)
        # PoolingCreateControl cbind converts the new vector to character.
        counts = np.array([float(format(x, ".15g")) for x in counts])
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.log2(test[:, 5].astype(float) / counts)
    if not np.isfinite(result).all():
        raise ValueError("legacy ratios contain non-finite values")
    classes = test[:, 6]
    if not np.isin(classes, ["IN", "OUT"]).all():
        raise ValueError("unrecognized target class")
    for kind in ["IN", "OUT"]:
        mask = classes == kind
        if mask.any():
            result[mask] -= np.median(result[mask])
    return result


def segment_profile(matrix, values, target, parameters):
    h, fc = parameters["HSLM"], parameters["FastCall"]
    estimated = estimate_parameters(values, float(h["Omega"]))
    rows, paths = [], []
    indices = []
    for chromosome in target["chromosomes"]:
        selected = np.flatnonzero(matrix[:, 0] == chromosome)
        if len(selected) < 2:
            raise ValueError(f"unsupported empty or single-window chromosome: {chromosome}")
        positions = np.trunc(matrix[selected, 1].astype(float))
        first, last = target["centromeres"][chromosome]
        before, after = np.flatnonzero(positions < first), np.flatnonzero(positions > last)
        arms = [np.arange(len(selected))]
        if len(before):
            if not len(after):
                raise ValueError(f"legacy arm split has no long arm: {chromosome}")
            arms = [np.arange(before[-1] + 1), np.arange(after[0], len(selected))]
            if sum(map(len, arms)) != len(selected):
                raise ValueError("legacy arm split excludes windows inside the centromere")
        for arm in arms:
            index = selected[arm]
            fit = segment(
                values[index],
                positions[arm],
                estimated,
                theta=float(h["Theta"]),
                distance=float(h["D_norm"]),
                min_windows=fc["minExons"],
            )
            for i, value in zip(index, fit.values, strict=True):
                rows.append(
                    [
                        matrix[i, 0],
                        str(int(float(matrix[i, 1]))),
                        number(matrix[i, 2]),
                        number(matrix[i, 3]),
                        number(values[i]),
                        number(value),
                        matrix[i, 6],
                    ]
                )
            paths.append(fit.path)
            indices.extend(index)
    if len(indices) != len(values) or len(set(indices)) != len(values):
        raise ValueError("target chromosome list does not cover every prepared window")
    return np.array(rows), np.concatenate(paths), np.array(indices)


def summarize(rows):
    """MakeData groups by consecutive segment value, even across chromosomes."""
    segmented = rows[:, 5].astype(float)
    starts = np.r_[0, np.flatnonzero(np.diff(segmented) != 0) + 1]
    ends = np.r_[starts[1:], len(rows)]
    return starts, ends, segmented[starts]


def run_analysis(
    samples,
    input_folder,
    target_folder,
    output,
    experiment,
    parameters=None,
    threads=1,
    force=False,
    *,
    r_seed_states=None,
):
    if threads != 1:
        raise ValueError("M3 supports --threads 1; parallel analysis is not yet qualified")
    if force:
        raise ValueError("M3 does not overwrite results; choose a new output directory")
    output = Path(output)
    if output.exists():
        raise ValueError("output directory already exists")
    input_folder, target_folder = Path(input_folder), Path(target_folder)
    prepared = load_manifest(input_folder, "prepared")
    target = load_manifest(target_folder, "target")
    if prepared["target_id"] != target["target_id"]:
        raise ValueError("prepared and target identities differ")
    design = experimental_design(read_yaml(samples), experiment)
    seeds = {}
    if r_seed_states is not None:
        seeds = read_yaml(r_seed_states)
        if set(seeds) != {test for test, _ in design}:
            raise ValueError("R RNG states must contain exactly the analyzed test sample names")
        # Validate every state before computing or publishing anything; empty
        # posteriors consume no RNG words and return the canonical signed state.
        seeds = {
            name: assign_labels(np.empty((0, 5)), r_seed=state).r_seed
            for name, state in seeds.items()
        }
    params = read_yaml(parameters) if parameters else DEFAULTS
    if set(params) != set(DEFAULTS) or any(set(params[k]) != set(DEFAULTS[k]) for k in DEFAULTS):
        raise ValueError("parameters must contain the documented HSLM and FastCall fields")
    matrices = {}
    for name in dict.fromkeys(name for test, controls in design for name in [test, *controls]):
        if name not in prepared["samples"]:
            raise ValueError(f"missing prepared sample: {name}")
        filename = prepared["samples"][name]
        if filename not in prepared["files"]:
            raise ValueError("sample file is not covered by manifest checksum")
        with np.load(input_folder / filename, allow_pickle=False) as archive:
            matrix = archive["matrix"]
        if matrix.ndim != 2 or matrix.shape[1] != 7:
            raise ValueError("expected seven-column prepared matrix")
        matrices[name] = matrix
    metadata = next(iter(matrices.values()))[:, [0, 1, 2, 3, 4, 6]]
    if any(not np.array_equal(m[:, [0, 1, 2, 3, 4, 6]], metadata) for m in matrices.values()):
        raise ValueError("prepared sample windows differ")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".excavator2-", dir=output.parent))
    try:
        summaries = {}
        for test, controls in design:
            values = ratios(matrices[test], [matrices[c] for c in controls])
            rows, path, indices = segment_profile(matrices[test], values, target, params)
            starts, ends, original = summarize(rows)
            fc = params["FastCall"]
            fit = fit_fastcall(
                correct_cellularity(original, float(fc["Cellularity"])),
                upper=float(fc["u"]),
                lower=float(fc["d"]),
            )
            calls = assign_labels(fit.posterior, r_seed=seeds.get(test))
            folder = temporary / "Results" / test
            folder.mkdir(parents=True)
            write_results(folder, test, rows, starts, ends, original, calls, target, target_folder)
            np.savez_compressed(
                folder / "checkpoints.npz",
                ratios=values,
                row_indices=indices,
                path=path,
                segment_starts=starts,
                segment_ends=ends,
                posterior=fit.posterior,
                labels=calls.labels,
                fastcall_trace=fit.trace,
                **(
                    {"r_seed_before": seeds[test], "r_seed_after": calls.r_seed}
                    if test in seeds
                    else {}
                ),
            )
            summaries[test] = {
                "windows": len(rows),
                "segments": len(starts),
                "calls": int(np.count_nonzero(calls.labels)),
                "fastcall_iterations": fit.iterations,
            }
        manifest = {
            "schema": 1,
            "kind": "analysis",
            "target_id": target["target_id"],
            "experiment": experiment,
            "parameters": params,
            "samples": summaries,
            "input_manifest_sha256": digest(input_folder / "manifest.json"),
            "target_manifest_sha256": digest(target_folder / "manifest.json"),
            "sample_design_sha256": digest(samples),
            "backend": "scalar",
            "threads": 1,
            "rng_policy": "per-sample-original-R-state" if seeds else "reject-random-ties",
            "r_seed_states_sha256": digest(r_seed_states) if seeds else None,
        }
        (temporary / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        temporary.rename(output)
    except BaseException:
        shutil.rmtree(temporary)
        raise
