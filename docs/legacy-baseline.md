# Legacy baseline established

On 2026-09-24, two independent runs of the original three-stage pipeline completed
on the repository's paired BAM fixture with one worker. All **476 compared files
matched exactly**: saved-stage RData exports and final HSLM/FastCall/VCF outputs.
Only the VCF date field was removed from comparison. Plots were generated but
were not compared byte-for-byte. The full runs did not inject an RNG seed.

The test sample produced **281,567 HSLM window rows and zero non-normal CNV calls**.
Both samples have counts for all 23 supported chromosomes and normalized outputs.
This is a useful repeatable baseline, but the absence of non-normal calls means
it cannot by itself validate every FastCall state. A separate 15-segment synthetic
fixture from the original R implementation exercises labels -2, -1, 0, 1, and 2.
It includes binary fitted parameters/posteriors and declares its fixed RNG seed.

## Reproducibility records

- [Pinned source/image and reference URLs](../tools/oracle/baseline.json).
- [Input checksums](../tools/oracle/inputs.sha256.json).
- [Baseline report](../tools/oracle/baseline-report.json), including run times and native-library hashes.
- [Golden checkpoint/result checksums](../tools/oracle/golden.sha256.json).
- [Exact installed Conda package URLs](../tools/oracle/environment-explicit.txt).
- [Five-state FastCall fixture](../tests/fixtures/legacy-fastcall/README.md).
- [Runner instructions](../tools/oracle/README.md).

Large files remain locally in `.oracle/run-02/` and `.oracle/run-03/`, with raw
outputs, manifests, stage logs, runtime inventories, and binary checkpoints.
They are ignored by Git. `.oracle/run-01/` preserves an initial failed compilation
attempt; it is not a golden run. Both successful runs use the frozen source commit,
not the current working-tree implementation. No legacy scientific bug was fixed.

## Runtime limitation and observed quirks

The pinned image's compiler is incomplete: rebuilding Fortran failed because its
compiler `specs` file is missing. The original bundled shared libraries load,
and their accompanying Fortran sources match the frozen repository exactly.
The successful runs use those image-pinned binaries, whose hashes are recorded.
This establishes a reproducible runtime baseline, not a reproducible binary build.

The preparation wrapper initializes MAPQ to 0 and its worker to 20, but the wrapper
does not forward that parameter. The literal awk predicate retained all tested
MAPQ values (0, 1, 10, 19, 20, 60). This behavior is characterized for compatibility
and remains unchanged. Details are in
[characterization.json](../tools/oracle/characterization.json).

The RData exporter was checked with binary64 fractions, signed zero, infinities,
separate NA/NaN masks, matrix shape, and character values. It avoids a decimal
roundtrip for numeric arrays. Unsupported object types cause export failure.

## Porting status

M2 is complete: the [FastCall implementation](fastcall-reference.md) passes the
saved synthetic and baseline-derived fixtures with both Python and scalar C++
kernels. Initialization, iteration/stopping, cellularity, and assignment remain
in Python. M3 has not started. Continue using the original pipeline for analyses
until end-to-end parity is demonstrated.

Pooled designs, CRAM, alternative thread counts, and per-iteration HSLM
instrumentation remain additional coverage. FastCall boundaries, ties and
per-iteration traces are covered by the M2 fixtures. This baseline does not yet
establish equivalence of an all-new pipeline.
