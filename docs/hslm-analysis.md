# M3: HSLM and analysis from legacy prepared data

M3 is implemented for single-profile paired and pooled analysis from converted
legacy normalized counts. Python owns the experimental design, ratios, global
parameter estimation, chromosome arms, filtering, segment summaries, FastCall
control and writers. C++ owns only HSLM transitions, emissions and Viterbi, plus
the existing FastCall kernels. [M4 now adds BAM preparation](preparation.md); target generation remains
unimplemented.

## Running analysis

First export the original prepared and target RData with `tools/oracle/export.R`
inside the pinned legacy environment. The existing baseline already has these
exports under `.oracle/run-03/checkpoints/`. Convert them once:

```sh
uv run python tools/convert_legacy.py \
  --prepared-exports .oracle/run-03/checkpoints/prepared \
  --target-exports .oracle/run-03/checkpoints/target/hg38/SureSelectV7/w_30000 \
  --chromosomes .oracle/run-03/target/hg38/SureSelectV7/w_30000/SureSelectV7_chromosome.txt \
  --centromeres .test/ref/CentromerePosition_hg38.txt \
  --assembly hg38 --output converted/

uv run excavator2 analyze --samples .oracle/run-03/design.yaml \
  --input converted/prepared/ --target converted/target/ \
  --output results/ --experiment paired
```

Analysis uses no R, Perl or Fortran runtime. Raw RData is not accepted directly.
`--experiment pooling` (or `pooled`) averages C-labelled controls in YAML order;
only T-labelled samples are analyzed, matching the original job splitter.
Paired analysis sorts numeric suffixes and requires matching C/T pairs. Sample
names must be plain path components; duplicate YAML keys/names are rejected.
`--parameters` accepts the original HSLM/FastCall YAML fields and defaults.
This initial implementation requires `--threads 1` and a new output directory;
`--force` is explicitly unsupported. Results are assembled in a temporary sibling
and published after every sample succeeds. Plots are not generated.

Outputs under `Results/<sample>/` retain the original HSLM/FastCall table and
region/window VCF names. `checkpoints.npz` also records ratios, row indices,
Viterbi paths, segment boundaries, posteriors, labels and FastCall traces.
`manifest.json` records input hashes, parameters, design and backend.

## Artifact boundary

Version 1 artifacts contain JSON manifests and compressed, non-pickled NumPy
arrays. Prepared matrices retain the original seven character columns, including
numeric strings; this avoids bypassing the original R character roundtrip.
FRB reference matrices, chromosome order and centromeres are captured in the
target artifact. File SHA-256 checks and a shared hash of ordered window metadata
prevent accidental mixing. Prepared metadata is compared again across samples.
The converter consumes only the explicit `export.R` format; it does not implement
a general R serialization reader. This is the minimal analysis interchange, not
a completed general artifact contract for future target/preparation stages.

## Following the algorithm

- `src/excavator2/hslm.py`: parameter estimation, state grid, distance covariates,
  breakpoints, short-segment filtering and median reconstruction.
- `src/excavator2/reference/hslm.py`: readable, literal O(T K²) numerical reference.
- `cpp/hslm.cpp`: the same scalar recurrence, with strict first-state tie handling.
- `cpp/bindings/hslm.cpp`: finite/aligned/contiguous binary64 validation and GIL release.
- `src/excavator2/analyze.py`: design, ratios, chromosome arms and calling orchestration.
- `src/excavator2/writers.py`: original table/VCF policy and rounding.

The public `segment(..., backend="python", trace=True)` path returns reference
transitions, emissions, scores and predecessors. Native mode uses rolling scores,
on-demand transition blocks and an O(T K) traceback buffer; `trace=True` exposes
the full matrices for comparison. Matrix orientation is `(time, source,
destination)` for transitions and `(time, state)` for other traces. Paths and
predecessors retain Fortran's one-based state numbers. Inputs are not mutated;
callers must not mutate shared arrays while a GIL-released kernel reads them.

Important preserved behaviors:

- Global 1–99% inclusive trimming with sample variance precedes arm splitting.
- The binary64 state grid follows R's multiplication/addition order.
- The original unsuffixed Fortran PI literal first rounds to binary32. Using a
  more precise PI shifted every emission by about 1.39e-8; the port reproduces
  the observed value instead.
- `ELNSUM` uses the original log/exp order, without substituting `log1p`.
- A sole short segment loses its final breakpoint and reconstructs as zero.
- `MakeData` groups consecutive equal segment values even across chromosomes.
- VCF CN rounds the two-decimal fractional CN; table CN rounds the unrounded CN.
- Reference bases retain the original FRB selection order. VCF source/date/header
  conventions are retained for compatibility, including the historical date order.

## Evidence and limits

Five small original-R/Fortran fixtures cover shifts, short segments, duplicate
positions, single-window numerical kernels and symmetric ties. Paths,
predecessors, filtered breakpoints and reconstructed values match exactly.
Transitions, emissions and scores use rtol/atol 2e-13; these tolerances never
permit changed discrete decisions. Single-window wrapper calls are rejected
because the original `SortState` wrapper fails on them.

The non-normal fixture exercises two samples, two chromosomes with both arms,
four calls per sample, paired/pooling designs and both VCF writers. All 16 output
files match the original byte-for-byte excluding VCF dates. Fixture generators,
source/image identities and hashes are retained under `tools/oracle/` and
`tests/fixtures/legacy-analysis/`.

On the supplied real paired dataset, all 281,567 HSLM rows, coordinates, segment
boundaries and call decisions agree. The maximum continuous-value difference is
about 1.01e-16. The FastCall table and both VCFs match exactly except the date.
This baseline has zero non-normal calls, which is why the synthetic fixture is
also mandatory. See `tools/oracle/m3-paired-report.json`; reproduce the comparison
with `tools/oracle/compare_analysis.py ORIGINAL REWRITTEN --report report.json`.

All 77 tests pass both in the editable install and against a rebuilt wheel
installed outside the checkout. Local validation uses Python 3.12.6 on macOS ARM64 with Apple Clang 17. The standalone
HSLM smoke test passes AddressSanitizer/UndefinedBehaviorSanitizer for 1, 17 and
10,003 windows. CI includes Linux/macOS Python 3.11–3.13 and installed-wheel tests;
remote qualification remains pending.

A 501-window, 21-state, three-repeat local benchmark measured 51.2 ms for the
literal Python reference and 0.257 ms for native (about 199x). It compares against
explicit Python loops, not the legacy Fortran or an optimized NumPy algorithm,
and is not a full-pipeline speedup claim. Run `python benchmarks/hslm.py`.
No SIMD/AVX or alternative recurrence is enabled.

Unsupported cases fail explicitly: constant/non-finite global variance,
single-window/empty arms, decreasing positions, centromere splits that exclude
windows or lack a long arm, and random FastCall ties without a recorded R state.
The lower-level FastCall API can replay an exported R RNG state; the CLI does not
yet accept one. Multi-profile HSLM, CRAM/preparation, plots, arbitrary malformed
legacy inputs and cross-platform numerical qualification remain outside this
M3 acceptance scope. No known scientific bug has been intentionally repaired.

## Live CI concordance

Every push and pull request now also runs a fresh pinned-legacy versus current-wheel
comparison for paired and pooled non-normal analysis. Reports, outputs and failure
logs are uploaded as the `legacy-concordance` artifact. See
[the oracle harness instructions](../tools/oracle/README.md#live-old-versus-current-ci-concordance)
for reproduction and the exact gates. The job now also covers original versus current preparation on small BAMs;
the large supplied BAM/reference run is not part of this per-PR job. The new harness passes locally; its first GitHub-hosted run
will occur when these commits are pushed.
