# M4: BAM preparation with legacy parity

> Historical parity record. Current preparation and analysis follow the
> [signal calibration contract](signal-calibration.md): preserve raw counts,
> require `--calibration`, and apply no covariate correction or median centering.
> Legacy prepared inputs and the analysis commands below are superseded.


M4 adds the `prepare` command against an exported original target. BAM decoding
uses pysam/HTSlib; batched counting uses NumPy searches, and size/MAP/GC correction
remains readable Python/NumPy. No custom C++ was added. Existing HSLM/FastCall kernels
continue to handle analysis. M5 now provides native target generation; see
[target usage and acceptance](target-porting.md#m54-target-cli-and-supplied-data-acceptance).

## Usage

Either generate a target with `excavator2 target`, or convert the original target once. The `export.R` checkpoint format must include
MyTarget, GCC, MAP and FRB; pre-existing normalized samples are no longer required.
For the saved supplied-data baseline:

```sh
uv run python tools/convert_legacy.py \
  --target-exports .oracle/run-03/checkpoints/target/hg38/SureSelectV7/w_30000 \
  --chromosomes .oracle/run-03/target/hg38/SureSelectV7/w_30000/SureSelectV7_chromosome.txt \
  --centromeres .test/ref/CentromerePosition_hg38.txt \
  --assembly hg38 --output converted/

uv run excavator2 prepare --samples .test/sample_sheet.yaml \
  --target converted/target/ --output prepared/ --threads 4

uv run excavator2 analyze --samples .oracle/run-03/design.yaml \
  --input prepared/ --target converted/target/ --output results/ --experiment paired
```

Preparation YAML maps sample names to indexed BAM paths. Relative paths follow
the original convention: resolve them from the working directory. Independent
samples can run concurrently; each chromosome's selected-position stream and
normalization operation order remain fixed. One- and four-thread fixture runs
produce identical prepared matrices. Analysis still requires one thread.

Outputs retain the version 1 prepared manifest accepted by `analyze`. Each sample
NPZ stores its original-style seven-column character matrix plus integer counts,
weighted counts and the size/MAP/GC/final normalization checkpoints. The manifest
records BAM hashes, target/sample-manifest hashes, requested thread/MAPQ settings,
chunk size and the actual filter policy. Partial work is removed on failure;
existing outputs and `--force` are rejected. Plots are not generated.

Targets converted during M3 remain usable for analysis. Preparation requires a
new conversion containing `preparation.npz` with target rows, GC and mappability.
The converter verifies the shared ordered-window identity when both original
prepared samples and target features are supplied.

## Current counting policy

Preparation excludes flags 4 and 1024 and requires MAPQ >= `--mapq` (default 20).
Other alignment flags retain their existing treatment. Positions are one-based
SAM starts, counted independently in each inclusive `[start,end]` interval.
Overlapping/nested windows can count the same read. Operational chunks contain
500,000 selected records, but partitioning does not change counts. Empty streams
produce zeros; terminal windows and exact multiples of the chunk size are handled.

These replace P01/P04/P05/P08 (#21–#23, #25). `reference/reads.py` is retained only
as a historical legacy oracle, not a supported counting mode. Regenerate prepared
artifacts to obtain corrected counts. Manifests record the filtering and counting
policies. Worker failures include the sample name and cause a nonzero CLI exit;
no partial preparation is published (P12, #27).

`normalization.py` preserves raw density (`counts / (end-start)`), including zeros.
Calibration belongs to the test/reference analysis stage. See
[target corrections](target-porting.md) for geometry and feature changes.

## Current verification

Tests compare counts with an independent interval-membership calculation across
legacy defect fixtures and arbitrary partitions, verify MAPQ filtering through the
CLI, and exercise atomic failure with both sequential and concurrent workers.

## Historical verification

- All 281,567 integer counts agree exactly for each supplied BAM (563,134 window
  comparisons). Final prepared matrices agree exactly, including numeric strings.
  See `tools/oracle/m4-preparation-report.json`.
- New `prepare → analyze` retains all original coordinates, segment decisions,
  FastCall calls and VCFs on the supplied paired dataset. Continuous HSLM drift
  remains at about 1.01e-16. See `tools/oracle/m4-paired-report.json`.
- Original R/Fortran fixtures cover endpoints, gaps, nested windows, final-only
  reads, and 499,999 / 500,000 / 500,001 / 1,000,001 selected-read cases. Correction
  fixtures cover zero medians, bin edges, out-of-range covariates and short final bins.
- Live CI now starts from tiny BAMs as well as the M3 prepared-data fixtures. It
  executes the pinned filtering/counting/normalization scripts, compares exact
  counts and five numerical checkpoints, and compares paired/pooling call outputs.
  This path passed locally. Reports and both output sets are uploaded on failure.
- All 98 tests pass in editable and installed-wheel environments, including
  threaded preparation and failed-run cleanup.
  Cross-platform qualification of the new preparation path awaits CI after pushing.

On the supplied 281,567-window Test1 BAM, a local single-thread measurement took
about 0.20 seconds for decoding/counting and 0.52 seconds for correction on macOS
ARM64/Python 3.12.6. These exclude artifact serialization and are not comparisons
against legacy runtime. Reproduce with `python benchmarks/preparation.py BAM TARGET`.
These timings do not justify another custom native kernel.

CRAM, unindexed/unsorted BAMs, empty selected streams and undefined non-finite
normalizations are unsupported and fail explicitly. Full reference/target
construction remains M5; the per-PR live job uses small synthetic inputs rather
than downloading the full supplied BAM/reference dataset.
