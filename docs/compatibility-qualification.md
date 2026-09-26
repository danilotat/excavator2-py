# M6 compatibility qualification

M6 qualifies the Python/C++ compatibility implementation for the supplied BAM
acceptance dataset and the committed differential corpus. The user explicitly
deferred plotting on 2026-09-26; scientific arrays, decisions, tables and VCFs are
the acceptance scope. No original scientific defect is repaired by this milestone.

M6.1 began at `ac53247279f28966437c7788c49caa0196064941`; its evidence below is
historical. The final qualification section distinguishes the final code's checks
from those earlier runs.

## Remote evidence

[Python port run 36102519608](https://github.com/danilotat/excavator2-py/actions/runs/36102519608)
passed all seven jobs for that exact commit:

- Linux (`ubuntu-latest`), Python 3.11, 3.12 and 3.13.
- macOS (`macos-latest`), Python 3.11, 3.12 and 3.13.
- Live original-versus-current concordance for geometry, reference features,
  BAM preparation, paired analysis and pooled analysis.

Each package job runs lint, formatting, fixture tests, native sanitizer smoke
checks, wheel builds and tests against the installed wheel outside the checkout.
The live job runs the digest-pinned original runtime and current wheel on fresh
small inputs. Success establishes the configured runner matrix; it does not imply
Windows, every CPU architecture, every dependency version or alternate SIMD paths
are qualified.

## Independent clean-wheel acceptance

A new Python 3.12 environment was populated from a freshly built wheel and its
runtime dependencies, with no editable installation. Package and extension paths
were verified under that environment's `site-packages`. Commands ran from outside
the checkout with Python `-I` isolation. Only input paths in the example YAML were
resolved to absolute paths; scientific settings and input files were unchanged.

The verification covers all 137 tests and a new supplied-data
`target → prepare → analyze` paired run. Target arrays compare directly with the
converted original target; preparation arrays compare with the M4 artifacts
already qualified against the original; final analysis compares directly with the
saved original results. Reports and runtime/build provenance are retained in
`tools/oracle/m6-clean-wheel-report.json`.

To reproduce, create a fresh environment, install the built wheel plus pytest,
resolve input paths for execution outside the checkout, and run:

```sh
/path/to/clean/bin/python -I -m pytest /absolute/repository/tests
/path/to/clean/bin/python -I -m excavator2 target --config config.yaml --output target/
/path/to/clean/bin/python -I -m excavator2 prepare --samples samples.yaml --target target/ --output prepared/ --threads 4
/path/to/clean/bin/python -I -m excavator2 analyze --samples design.yaml --target target/ --input prepared/ --output results/ --experiment paired
```

Use `tools/oracle/compare_target.py` and `compare_analysis.py` for the exact target
and decision-gated analysis comparisons. Compare every preparation NPZ array,
not just final calls. See the M5 reports for the original comparison basis.

## Qualification matrix and remaining work

| Area | Evidence established | Remaining gap |
| --- | --- | --- |
| Supplied paired pipeline | New target through final calls; 281,567 windows per sample | Supplied example has no non-normal calls |
| Non-normal paired/pooling | Original fixtures, stage CI and M6.2 complete fresh-reference chain | Remote confirmation of the newly added full-chain CI step |
| Target geometry | Eight cases; live old/new comparison; supplied target exact | Additional assemblies, unusual contigs and broader interval edge cases |
| Target features | Ten cases plus all supplied arrays; binary32/R decimal replay; exact parity immediately below/at UCSC's 3,000-block switch | Additional assemblies and broader missing-data patterns |
| Read counting | Endpoints, overlaps, terminal/chunk cases; supplied exact counts | CRAM explicitly unsupported; empty streams rejected; unmapped/duplicate filtering and secondary/supplementary/QC-failed retention tested |
| Numerical kernels | Original fixtures, native/reference comparison, sanitizer smoke checks; M6.4 pure-R boundaries; M6.5 full-inference centromere cases | Original near-tie replay and empty-chromosome rejection covered; arbitrary invalid legacy domains excluded |
| Reproducibility | One/four-worker preparation tests; pinned oracle; installed wheels | Per-sample CLI RNG replay now covered; analysis parallelism remains unsupported |
| Plots and diagnostics | Scientific data tables/VCF comparisons; [original output inventory and acceptance gates](plotting-scope.md) | Deferred explicitly by user for M6; no visual parity claim |
| Packaging | Linux/macOS Python 3.11–3.13 CI and independent clean-wheel run | Release notices, migration guide and distribution qualification belong to M8 |

## M6.3 reference-feature boundary qualification

M6.3 qualifies the UCSC 3,000-block algorithm switch with generated 2,999- and
3,000-row inputs. The pinned original log proves that the cases exercise the
per-interval and chromosome-buffer paths respectively. Both use the same cycle of
fully covered, partially covered, zero-valued and uncovered intervals; their 2,999
shared MAP values are identical, and the port exactly matches every GC, MAP and
FRB value on both sides of the switch.

M6.4 adds the bounded numerical edge coverage below. Remaining numerical cases
and plot scope still need resolution before signing off M6. Keep scientific bug
fixes and performance changes separate from qualification.

## M6.4 HSLM wrapper boundary qualification

Fresh execution of the pinned original `LibraryJSLMIn.R` establishes nine
segment-filtering cases: sole segments below/at/above the four-window threshold,
short first/middle/last segments, adjacent short segments, all-short segments, and
unfiltered segments. Filtered breakpoints and reconstructed medians match exactly.
The original sole-short-segment zero reconstruction remains preserved.

Five parameter-estimation probes distinguish a valid varying profile from
constant, singleton, two-value and three-value profiles. Original percentile
trimming yields `NA` deviations when fewer than two observations remain, and zero
deviations for a constant profile. The port explicitly rejects those unsupported
domains rather than passing invalid parameters into native code. Original
`SortState` also fails on a single window; both port backends reject that wrapper
input even when global parameters are valid. These probes exercise pure R only;
they do not claim defined Fortran behavior for invalid parameters or qualify all
centromere/arm geometries.

Captured TSVs, the original error, source and harness hashes, and pinned runtime
identity are in `tests/fixtures/legacy-hslm/edges/`. Regenerate and compare with:

```sh
.venv/bin/python tools/oracle/characterize_hslm_edges.py \
  --output .oracle/hslm-edges-check --compare-fixtures
```

The live concordance CI job now runs this check; fixture-based Python checks run
in every package job. No production algorithm or original source was changed.
Local verification passed: 163 tests, lint, formatting and diff checks, plus an
independent fresh Docker replay with byte-identical boundary outputs. The new CI
step has not yet been run remotely.

## Separation of legacy issues and port evidence

The local, intentionally uncommitted `LEGACY_ISSUES.local.md` records original
bugs, hardcoded parameters, limitations and questionable assumptions only. Port
defects, qualification results and oracle-environment limitations belong in these
tracked documents. In particular, bedtools binary32 GC rounding and the pinned
R decimal-conversion replay are compatibility requirements whose port fixes are
documented in `target-porting.md`; the UCSC 3,000-block switch is qualified above,
and the original container's incomplete compiler is documented in
`legacy-baseline.md`. Mismatched original feature rows remain a legacy issue;
the port's artifact rejection is a validation policy, not another original bug.

## M6.5 chromosome-arm boundary qualification

Eight centromere configurations run through the full pinned original inference
script on the existing non-normal two-chromosome prepared fixture. Three succeed:
both arms, no short arm, and no short arm with some windows inside the centromere.
For both test samples, the port matches the resulting HSLM metadata, window
retention, segment boundaries and continuous values under the existing tolerances.
The asymmetric no-short-arm behavior remains preserved.

Five configurations fail in original inference: no long arm, a single-window
short arm, a single-window long arm, windows exactly on centromere endpoints,
and windows inside a two-arm centromere split. Original logs and exit statuses
are captured; the port rejects these geometries explicitly. Singleton cases
place the centromere between observed windows to isolate that failure from
centromere-window exclusion. Original failed runs stop at Test1; port rejection
is checked for each prepared test sample.

Evidence and hashes are in `tests/fixtures/legacy-analysis/arms/`. Reproduce with:

```sh
.venv/bin/python tools/oracle/characterize_arms.py \
  --output .oracle/arms-check --compare-fixtures
```

CI regenerates this evidence in its live concordance job. This step qualifies
arm splitting and HSLM output, not additional final-call or plot parity.
Local verification passed: 180 tests, lint, formatting and diff checks, and an
independent Docker replay with identical evidence files. Remote execution of
the new CI step remains pending.


## M6.2 fresh-reference non-normal chain

The installed-wheel run now covers the complete original and current target,
preparation and analysis CLIs on generated inputs: 23 reference contigs, 4,113
windows, two controls and one test with decreases/increases on chr1 and chr2.
Both paired and pooling produce four non-normal calls (losses and gains).

All target/GC/MAP/FRB arrays, all 12,339 raw counts and all prepared matrices match
exactly. Each mode's four scientific output files matches exactly apart from VCF
fileDate; this fixture has zero measured continuous HSLM drift. The harness rejects
missing output files, fewer/different calls, and loss of either loss/gain coverage.
It requires exact target and preparation equality before accepting final calls.

See `tools/oracle/m6-end-to-end-report.json` and reproduce using:

```sh
/path/to/installed-wheel/python tools/oracle/end_to_end.py --output .oracle/end-to-end
```

The current path generates its own target; converted original artifacts are used
only as comparison inputs. No downloaded genome/BAM assets or previous oracle
runs are required. CI now invokes the same harness and retains its complete logs,
inputs and report. Local acceptance passed; the new remote CI step is pending.
Six comparator tests additionally prove that valid-checksum scientific changes
in target features, row ordering, reference bases and centromeres are rejected.
The full local suite has 143 passing tests.

Plot qualification remains open: the original plotting code cannot obtain UCSC
chromInfo for the synthetic assembly in the offline container and fails in both
modes. This is explicit in the report; the pass applies to scientific data and
calls only. A first fixture with longer overlapping IN windows also exposed
nonmonotonic midpoints and failed in both original inference and port validation.
The passing fixture uses shorter intervals, without changing either scientific
implementation. Reproducer details are retained in the local improvement log.

## Final M6 acceptance scope

The supported production path uses indexed BAM, the legacy chromosome/reference
policy, scalar binary64 kernels, one analysis worker, and either unique FastCall
maxima or recorded per-sample original R states. Preparation is qualified at one
and four workers. Successful original inputs in the committed corpus must retain
exact decisions; known undefined/degenerate domains fail explicitly. CRAM,
arbitrary assemblies/contigs, additional platforms, analysis parallelism and
alternative SIMD kernels are not silently included in the compatibility claim.
They are extensions requiring their own evidence, not unexplained mismatches.

The final CLI addition accepts `--r-seed-states` with original per-sample
Mersenne-Twister states. Tests replay recorded exact/near ties through analysis
orchestration, check independent sample streams, preserve unique-call outputs,
validate state bounds and provenance, and reject incomplete mappings before
publishing results. See `hslm-analysis.md` for the input contract. This is replay
of the original random decision rule, not a new tie-breaking policy.

Closure requires the full fixture suite in editable and clean-wheel environments,
the live old/new full chain, native sanitizer smoke checks, and successful remote
Linux/macOS Python 3.11–3.13 package/concordance jobs. Final evidence is recorded
in `tools/oracle/m6-final-report.json`. Plotting remains deferred as requested;
performance and release/distribution work remain M7 and M8 respectively.

### Final local result — 2026-09-26

Revision `06ea804` passes all 188 tests in both editable and isolated clean-wheel
runs, lint/format checks, and native AddressSanitizer/UndefinedBehaviorSanitizer
smoke tests. The installed Python sources match the committed production files.
A fresh supplied-data target → prepare → analyze run matches all 281,567 windows
and every prepared array against the qualified baseline; scientific outputs
match the saved original under the existing field-specific gates. A fresh
pinned-original synthetic run matches all target/preparation arrays and four
non-normal calls in each of paired and pooling. The M6.4/M6.5 boundary evidence
also reproduces independently. See `tools/oracle/m6-final-report.json` for artifact
hashes, runtime, source identity and complete comparison results.

**Local M6 acceptance is complete with plotting deferred. Final remote sign-off
is pending.** The previous revision `b08c9db` passed
[Python port CI](https://github.com/danilotat/excavator2-py/actions/runs/36225567930).
It does not qualify the new code. Automatic approval review blocked the push to
`danilotat/excavator2-py` because authorization for that destination was not
established; permission has been requested. Once authorized, push `dev/porting`,
require all package/concordance jobs to pass for the final code revision, and
record that run before marking the remote gate complete.
