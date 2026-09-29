# M4: BAM preparation with legacy parity

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

## Compatibility decisions

`reads.py` preserves the source's behavior, including its known defects:

- Only flags 4 and 1024 are excluded. Secondary/supplementary alignments, mate
  records and QC-failed records are retained. Coordinates use one-based SAM POS.
- The original awk expression retains every MAPQ. The CLI accepts `--mapq` and
  records it, but does not silently introduce a quality filter in compatibility mode.
- Chunks contain 500,000 selected records. A wholly contained chunk replaces the
  residual count instead of adding to it. An exact full final chunk can lead to
  the original empty-read failure; the port reports an error in that case too.
- A window is committed only when a later read advances beyond its end. An
  unflushed final window stays zero. Starts/ends are inclusive for counting.
- Overlapping/nested target windows are not repaired. Running maximum ends and
  search bounds preserve the original forward-only assignment. A literal Python
  state-machine reference is in `reference/reads.py` for debugging.
- The original Fortran accesses beyond its arrays after the last window. The port
  stops once all defined counts are written. This avoids unsafe memory access;
  all defined counts agree on the supplied data and fixtures.

`normalization.py` divides by `end-start`, then applies size correction to IN
windows, MAP correction separately to IN/OUT, and finally GC correction separately.
Bins are `[0,5]`, `(5,10]`, etc. A partial last size bin is omitted, matching R's
`seq`. Non-positive bin medians leave values unchanged. Final zeros are replaced
with the smallest nonzero value within the same IN/OUT class. The saved matrix
uses the original 15-significant-digit character conversion with `scipen=20`.
No pseudocount, alternative normalization, quality policy or scientific bug fix
is introduced.

## Verification

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
