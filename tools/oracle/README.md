# Pinned legacy oracle

`baseline.json` pins commit `92f9c6e79fdeae32a09b99fa6994edc9a53e4efc` and the
Linux/amd64 container image by digest. The runner extracts that commit's source,
not the current working tree or the image's Perl/R scripts. `inputs.sha256.json`
records the exact fixture/reference contents used for this baseline; new runs
reject changed inputs instead of silently replacing the baseline.

The pinned image has an incomplete Fortran compiler installation. Its original
precompiled libraries load successfully; the runner verifies that their
accompanying Fortran sources match the snapshot byte-for-byte, then uses those
binaries and records their hashes. This is a reproducible **runtime** baseline,
not a demonstrated rebuild of the legacy toolchain. No legacy scientific code
is patched.

## Prerequisites

- The development environment from `docs/development.md`.
- Docker running, with Linux/amd64 support (emulation on Apple Silicon).
- Download the FASTA gzip and BigWig from `reference_urls` in `baseline.json`
  into the paths in `.test/config.yaml`; decompress the FASTA.
- Generate the FASTA index once with `pysam.faidx` or `samtools faidx`.
- Pull the digest-pinned image in `baseline.json` before running. Execution itself
  uses `--network none` and mounts the repository read-only.

The large downloaded inputs and local `.oracle/` runs are ignored by Git. Keep
approximately 6 GB available for the downloaded/compressed/decompressed reference
assets, plus the container image and run outputs. Current `.test` inputs are BAM
files; this runner is intentionally scoped to that paired fixture.

## Run and compare

```sh
uv run python tools/oracle/oracle.py run --output .oracle/reference-1 --threads 1
uv run python tools/oracle/oracle.py run --output .oracle/reference-2 --threads 1
uv run python tools/oracle/oracle.py compare .oracle/reference-1 .oracle/reference-2 \
  --report .oracle/repeatability.json
```

Output directories must be new: failed runs are preserved for diagnosis and are
never overwritten. Each run records input/source/harness hashes, command arguments,
thread count, status, exit code, runtime package/tool information, and stage logs.
Future runs snapshot the harness into the output as well. The legacy parent scripts
can swallow worker failures, so the runner requires all expected per-sample NRC and
HSLM/FastCall/VCF result files before marking a run complete.

`export.R` exports all saved atomic RData objects after the unmodified pipeline:

- Double values: binary64, little endian, R column-major order.
- Integer/logical values: signed 32-bit little endian.
- Character values: JSON strings, retaining R's already-performed conversions.
- Separate missing/NaN masks and shape/dimname metadata.

Unsupported R object types fail export rather than being silently coerced.
These are saved-stage checkpoints; per-iteration HSLM/EM instrumentation remains
future work, and must be verified against the untouched outputs.

The comparator requires complete runs and matching source/input/runtime baseline
identity. It compares exported arrays and final text/VCF records exactly, preserving
row order. Only VCF `fileDate` is stripped. Plots, logs, and serialized RData bytes
are excluded from equality checks. This deliberately strict repeatability check
is not yet the field-tolerant comparator for a future Python/C++ implementation.
It reports every differing artifact and exits nonzero on any difference.

No RNG seed is injected into the original run. If repeated calls differ, investigate
and capture RNG state before accepting goldens. Run another repetition at `--threads 3`
after same-thread repeatability succeeds. A completed baseline is not evidence that
the port produces equivalent calls: no scientific algorithms have been ported yet.

## Additional characterization

`characterization.json` records the observed MAPQ predicate and the exporter smoke
check. `characterize_fastcall.R` generates a small five-state fixture from the original
R implementation with a declared fixed seed; its binary arrays and provenance are
stored in `tests/fixtures/legacy-fastcall/`. This complements the shallow paired
example, whose first reference run yielded no non-normal calls.

## Recorded result

Two full one-worker paired runs matched across 476 exported checkpoint/result files,
with only VCF dates excluded. See `baseline-report.json`, `golden.sha256.json`, and
`docs/legacy-baseline.md`. `environment-explicit.txt` records the image's installed
package URLs. Other experimental modes and worker counts remain unqualified.

## Live old-versus-current CI concordance

The `Live legacy concordance (paired and pooling)` job in `python-port.yml`
runs on every push, pull request, and manual workflow dispatch. It builds the
current wheel and tests it in a separate environment. Each run:

1. Extracts the original commit in `baseline.json` with `git archive`.
2. Runs that snapshot in the digest-pinned legacy container, with networking
   disabled inside the container. The bundled Fortran binary is used only after
   its accompanying source matches the snapshot. The container's R wrapper is
   not substituted for the pinned source.
3. Generates fresh two-chromosome normalized-count inputs, executes the original
   paired and pooled analyses, and exports those same inputs for the new package.
4. Runs the current wheel's CLI and compares all eight output files per design.

Calls and VCF records must match exactly (except VCF `fileDate`). HSLM coordinates,
classes, row order and segment boundaries must agree exactly; continuous HSLM
values retain rtol=1e-13 and atol=2e-15. Missing or extra output files fail the job.
Failures return a nonzero exit code; no `continue-on-error` or conditional skip
is used for concordance. The `legacy-concordance` artifact retains JSON reports,
source/input provenance, both sets of outputs and logs for 14 days, including
failed runs. The existing cross-platform fixture tests remain separate.

Reproduce with Docker available and the current package installed:

```sh
python tools/oracle/concordance.py --output .oracle/live-concordance
```

Choose a new output directory for each run. This performs real old-versus-new
execution, not just a comparison against committed golden files. It now covers M3 analysis and M4 preparation: small indexed BAMs are generated,
the original filtering/counting/normalization scripts run, and exact counts plus
five correction checkpoints are compared before paired/pooling calls. This adds
another 16 call-output comparisons. It does not download large BAM/reference
assets or rerun target generation; the full supplied-data baseline remains separate.

## M5.1 target geometry characterization

Run the unchanged baseline `FilterTarget.R` on eight small synthetic cases:

```sh
python tools/oracle/characterize_target.py --output .oracle/target-geometry
```

Use a fresh output directory and the installed port environment with Docker.
Five cases succeed and three intentionally reproduce original failures. The
harness checks exported RData against `Filtered.txt`, then writes compact arrays,
inputs, stderr and provenance under `fixtures/` in the output directory. Committed
captures live in `tests/fixtures/legacy-target/`; their tests run in wheel CI.
The Python geometry port now compares every output cell against these captures
in the fixture tests. FASTA/BigWig features and live old-versus-new target CI are
later M5 steps. See
[the target contract](../../docs/target-porting.md).
