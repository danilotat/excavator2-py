# EXCAVATOR2 Python 3 / C++ porting roadmap

Planning baseline: repository commit `92f9c6e79fdeae32a09b99fa6994edc9a53e4efc`, inspected on 2026-09-24.

This is an implementation plan, not a claim that the port or equivalence tests already exist. The legacy pipeline was inspected but not executed for this document.

## 1. Objective and release contract

Build an installable, Python-first package with a small C++ acceleration core. Readability and maintainability for a Python developer are architectural requirements. Python owns the pipeline, target construction, file access through existing libraries, normalization, experimental design, algorithm setup/control, segment filtering, calling policy, outputs, and plots. Custom C++ is reserved for demanding numerical kernels, initially HSLM transitions/emissions/Viterbi and FastCall's repeated E/M calculations.

Use Python with NumPy/SciPy and existing native-backed I/O libraries before writing additional C++. Array operations already execute compiled routines; a Python module containing numerical work does not imply slow Python element-by-element execution. Add a custom native kernel only when profiling identifies a material bottleneck or differential tests establish a numerical compatibility need. Keep the scientific choices visible in Python.

The first release must reproduce the original calls on the acceptance dataset before any scientific or behavioral bug is corrected. Preserve observable quirks in a named compatibility implementation; record suspected bugs separately. A mathematically cleaner implementation that changes a breakpoint or call fails this milestone.

Use the repository's `.test` paired dataset as the initial acceptance fixture. The user's “given file” is not otherwise identified here: if it refers to another BAM/CRAM or existing result set, add that exact input and its supporting files to the mandatory acceptance manifest before declaring parity. Do not substitute a synthetic dataset for it.

Define compatibility at three levels:

1. **Exact decisions:** target rows and order, integer read counts, retained/missing windows, HSLM paths and breakpoints, filtered segments, called intervals, discrete labels, rounded copy numbers, and VCF record identities must match.
2. **Numerical agreement:** continuous intermediate values and probabilities must meet explicit, field-specific tolerances established from reference runs. Such tolerances never excuse a changed decision.
3. **Artifact agreement:** preserve result filenames, columns, record order, coordinate conventions, missing-value spelling, and meaningful VCF fields. Timestamps and run locations may differ under a narrow documented whitelist. Internal RData byte compatibility and pixel-identical plots are not required.

The normal installed package must eventually run all three commands without R, Perl, Fortran, or shell pipelines. Those dependencies remain in an isolated development oracle and optional conversion tools. Temporary mixed-language implementations are development milestones, not the finished port.

Explicitly defer changed calling models, corrected interval semantics, added chromosomes, revised filters, float32 arithmetic, GPU execution, and new scientific defaults until the compatibility release has passed.

## 2. What the current implementation actually does

Paths in this document are relative to the repository root. The executable code, followed by recorded behavior of a pinned reference build, is the specification when README prose differs.

| Stage | Current implementation | Port responsibility |
|---|---|---|
| Target command | `excavator2/TargetPerla.pl` | Python settings, provenance, output orchestration |
| Target geometry | `lib/R/FilterTarget.R` beneath `excavator2/` | Python/NumPy interval construction, sorting, gap filtering |
| Reference features | `lib/bash/TargetCreate.sh`, `lib/R/SaveMap.R`, `SaveGCC.R`, `SaveFRB.R` | Python using existing FASTA/BigWig readers and array operations |
| Prepare command | `EXCAVATORDataPrepare.pl`, `lib/perl/ReadPerla.pl`, parallel job scripts | Python sample scheduling |
| Alignment selection and counts | `lib/bash/FiltBam.sh`, `lib/R/MakeReadCount.R`, `lib/F77/F4R.f:EXOMECOUNT` | Python reader/filter policy and chunk orchestration; small C++ counting kernel only if measured necessary |
| Normalization | `lib/R/EXCAVATORNormalizationExome.R`, `LibraryExomeRC.R` | Python/NumPy binning, medians, corrections, diagnostics |
| Controls and ratios | `PoolingCreateControl.R`, `EXCAVATORInferenceExome.R` | Python/NumPy pooling, transforms, and design mapping |
| HSLM | `LibraryJSLMIn.R`, `lib/F77/FastJointSLMLibraryI.f` | Python parameters, arm splitting, filtering; C++ transitions/emissions/Viterbi |
| FastCall | `LibraryFastCall.R`, `EXCAVATORInferenceExome.R` | Python initialization, EM loop, stopping/assignment policy; C++ E/M and posterior kernels |
| Results and plots | VCF routines in `LibraryFastCall.R`, plotting scripts | Python compatible serialization, summaries, and plotting |

The following observations deserve explicit characterization tests before translation:

* `FilterTarget.R` constructs chromosomes 1–22 and X explicitly, uses a 200-base flank, and applies its own gap test. Do not replace that test with a conventional overlap predicate without proving equivalence. Input ordering and chromosome naming are significant.
* Target feature extraction uses different coordinate transformations for GC/mappability and reference bases. `SaveMap.R` selects column six from `bigWigAverageOverBed`. Capture the tool's actual covered/missing-base behavior rather than assuming a generic interval mean.
* `FiltBam.sh` uses `samtools view -F 1028` and selects alignment start positions. This is not a base-coverage or fragment-overlap counter. Its awk expression contains a literal `"/$MAPQ/"`; characterize its effective comparison in the pinned shell/awk environment instead of implementing the apparent intended MAPQ rule.
* Additional oracle setup observation: the preparation wrapper initializes MAPQ to 0, while its worker initializes it to 20; the parallel-job generator does not forward the wrapper setting. A pinned-runtime probe retained MAPQ 0, 1, 10, 19, 20, and 60 through the literal awk predicate. Preserve this observed behavior until post-parity bug fixes.
* Counting reads proceeds in chunks of 500,000, carries residual state, and has potentially problematic terminal behavior. Counts may depend on wrapper behavior as well as the Fortran routine.
* Normalization divides by `end - start`. IN windows undergo size, mappability, and GC correction; OUT windows undergo mappability and GC correction. Corrections use bin medians, with step 5 and specific endpoint inclusion rules.
* Zero normalized values are replaced separately in IN and OUT groups with their smallest nonzero values. Empty and all-zero cases require characterization.
* Pooling code averages already normalized control values, despite the README describing summation. Preserve the actual operation and addition order.
* The Perl CLI accepts `pooling`, whereas parts of the README say `pooled`. The new CLI can accept both spellings and map them to the same legacy behavior.
* Log ratios are centered by separate IN/OUT medians. A comment mentions LOWESS, but the inspected inference path performs median centering; do not add smoothing.
* HSLM estimates variance after inclusive 1st/99th percentile filtering, then estimates parameters once before chromosome/arm segmentation. It uses 21 states for the single-profile path, from -1 to 1 in increments of 0.1.
* `SegResults` assigns segment **medians**, despite the output column being called `SegMean`.
* `FilterSeg` acts on segment window counts with `<= minExons`, and removes breakpoints using particular first-segment logic. It is not simply “drop segments containing fewer than four IN exons.”
* `MakeData` identifies segments through changes in successive segment values. Explicit chromosome boundaries are not part of that expression; equal adjacent values across boundaries need a regression fixture.
* FastCall keeps component means fixed; only deviations and priors are updated. Initialization ranges and E-step truncation ranges differ. Its stopping statistic is not a conventional mixture log-likelihood.
* `LabelAss` calls R's `max.col` without a tie-method argument. Characterize ties, near ties, and RNG dependence. A deterministic first-maximum replacement is not automatically equivalent.
* `D_norm: 10e5` in the supplied YAML is 1,000,000. Some README prose says 100,000. Ship the executable configuration value as the compatibility default.
* VCF copy-number rounding differs from the text path: VCF rounds fractional CN before rounding integer CN. Preserve these separate serialization paths.

These are observations and investigation targets, not a declaration that every item is a confirmed bug.

## 3. CLI and configuration

The primary interface is exactly the requested three-stage workflow:

```sh
# Initialize target regions
excavator2 target --config config.yaml --output target/

# Prepare read counts
excavator2 prepare --samples samples.yaml --target target/ --output prepared/ --threads 4

# Call CNVs
excavator2 analyze --samples sample_list.yaml --input prepared/ --target target/ --output results/ --experiment paired
```

| Command | Required options | Additional first-release options |
|---|---|---|
| `target` | `--config`, `--output` | `--force`, `--verbose` |
| `prepare` | `--samples`, `--target`, `--output` | `--threads 1`, `--mapq 20`, `--force`, `--verbose` |
| `analyze` | `--samples`, `--input`, `--target`, `--output`, `--experiment` | `--threads 1`, `--parameters`, `--force`, `--verbose` |

Provide `--help` and `--version`. Retain sensible original short aliases: `-s`, `-o`, `-t`, `-i`, `-e`, `-p`, `-q`, `-@`, `-v`, and `-f`, only where relevant. `target -s` aliases `--config`; accepting `--settings` eases migration. Legacy script-name launchers are optional and should only forward arguments to the new CLI.

Keep the three existing YAML schemas:

* Target settings: `Reference.{Assembly,FASTA,BigWig,Chromosomes,Centromeres,Gaps}` and `Target.{Name,BED,Window}`.
* Preparation: sample label to BAM/CRAM path.
* Analysis: experimental tag (`T1`, `C1`, etc.) to prepared sample label.

Keep the existing HSLM and FastCall parameter keys and defaults. Record the resolved values and parameter-file checksum in every analysis manifest. Normalize scientific-notation parsing explicitly.

Resolve relative paths against the invocation working directory for the initial CLI, matching the supplied `.test` examples; write resolved paths to provenance. Do not silently change to YAML-relative resolution. Handle paths with spaces without constructing shell command strings.

For `--target target/`, write a manifest directly at that root, even if legacy-compatible ancillary artifacts are nested beneath it. Resolve legacy target leaf directories through an explicit importer. Reject ambiguous directory discovery instead of selecting the first matching file.

Validate file existence, duplicate YAML keys, duplicate sample labels, missing prepared samples, reference identity, target identity, and paired design completeness. Preserve the reference pairing/order on valid legacy inputs. Distinguish unsupported malformed inputs from successful-but-quirky inputs that must retain compatibility.

Use stderr for progress and diagnostics. A nonzero exit must mean failure; partial output must not look complete. Use a per-output lock, run-specific temporary directory, and atomic manifest finalization. `--force` may replace owned artifacts in the selected output, never arbitrary external paths. Thread count is a global compute budget, including native decoding workers.

## 4. Package and native architecture

### 4.1 Proposed stack

Use C++17 for the initial native implementation, CPython with an explicitly tested version matrix, `pybind11` for bindings, and CMake through `scikit-build-core` for wheel/source builds. Select and pin concrete dependency versions during the packaging milestone, rather than depending on moving development versions.

`pybind11` supports typed NumPy arrays and buffer interfaces; these provide the proposed coarse-grained native boundary. Validate contiguity, dtype, shape, ownership, and strides deliberately rather than permitting unexpected copies. See the [official array/buffer documentation](https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html).

`scikit-build-core` is the proposed Python/CMake build backend. Keep build configuration in `pyproject.toml` and `CMakeLists.txt`, with reproducible release settings. See its [official documentation](https://scikit-build-core.readthedocs.io/en/latest/).

Use existing Python-accessible readers, with [pysam](https://pysam.readthedocs.io/en/latest/api.html) and [pyBigWig](https://github.com/deeptools/pyBigWig) as the initial candidates, for BAM/CRAM/FASTA and BigWig access. Validate their selected operations against the legacy tools before pinning versions. Keep EXCAVATOR2 filtering and coordinate policy in Python. Avoid maintaining a second custom C++ I/O layer. Batch processing and array operations should avoid unnecessary per-base work; measure per-alignment Python overhead before deciding whether a narrowly scoped extraction/counting accelerator is needed.

Use NumPy for arrays and vectorized calculations, SciPy where its numerical primitives pass compatibility tests, YAML parsing for configuration, and an optional plotting dependency. A large dataframe framework, general HMM package, or generic Gaussian-mixture implementation is unnecessary for the core and risks replacing legacy semantics.

### 4.2 Proposed layout

```text
pyproject.toml
CMakeLists.txt
src/excavator2/
    __init__.py
    cli.py
    config.py
    pipeline.py
    target.py                  # geometry and reference features
    reads.py                   # reader adapters, selection policy, chunks
    normalization.py           # readable NumPy correction steps
    design.py                  # pairing, pooling, log ratios
    hslm.py                    # parameters, arm splits, native call, filtering
    fastcall.py                # initialization, EM loop, stopping, labels
    compat.py                  # R conventions and conversion behavior
    reference/                 # small Python/NumPy versions of native kernels
        hslm.py
        fastcall.py
    artifacts.py
    provenance.py
    reporting.py
    _core.<extension>
cpp/
    include/excavator2/
    hslm.cpp                   # transitions, emissions, Viterbi and traceback
    fastcall.cpp               # repeated density/E/M/posterior calculations
    numerical_compat.cpp       # only helpers needed by these kernels
    simd/                      # optional; introduce after profiling and parity
    bindings/
tests/
    unit/
    differential/
    integration/
    fixtures/
    native/
tools/oracle/                  # pinned legacy runner and RData/checkpoint export
benchmarks/
docs/
```

Keep the original `excavator2/` tree intact as the inspected reference until the oracle is archived independently. The `src/` layout prevents installed Python code from being confused with that directory.

### 4.3 Boundaries and data ownership

Keep public stage APIs in Python. The initial private extension should expose only a few batch operations, such as `hslm_path`, `fastcall_estep`, `fastcall_mstep`, and `fastcall_posterior`. Python computes HSLM parameters and distances, passes arrays/scalars to the core, and reconstructs/filters segments. Python owns FastCall's initialization, iteration loop, stopping test, cellularity correction, labels, and rounding. Each FastCall kernel handles all segments in one call; never cross the binding once per observation. These names are proposals, not implemented functions.

Keep a straightforward Python/NumPy reference implementation of each custom kernel, usable on small fixtures through the Python API. It is both executable documentation and a debugging backend, not a promise of whole-genome Python performance. Document array shapes, equations, legacy function mapping, and one tiny worked example next to each kernel. Avoid template-heavy frameworks and general-purpose C++ class hierarchies.

Adding another C++ function requires a short benchmark or compatibility justification, an understandable Python reference, and a differential test. Keep filesystem paths, YAML, sample objects, output formatting, and scientific policy out of the extension. A future fused EM loop is optional only if profiling shows Python iteration overhead is material; preserve the Python loop as the specification.

Pass contiguous typed arrays, with explicit window IDs and group offsets. Return typed arrays and small diagnostics objects. Release the Python GIL during long native work; never access Python objects while it is released. Keep input arrays alive for the full call, and avoid hidden mutable global buffers. Raise translated exceptions with sample/chromosome context.

Use binary64 for numerical calculations. Use explicit-width coordinate/count types, preferably 64-bit for storage, while testing legacy 32-bit conversion boundaries. Do not rely on C++ signed overflow or reproduce memory corruption. Inputs whose reference results depend on undefined memory behavior must be documented as outside proven parity until resolved.

Keep C++ numerical code independently testable without Python. Put compatibility semantics in named functions with links to their legacy source and a regression test, so future scientific changes can be introduced deliberately.

### 4.4 Artifact schema

Use a versioned JSON manifest plus typed `.npy` arrays and small TSV metadata tables initially. Arrays should be loadable without pickle and memory-mappable where beneficial. Avoid one archive requiring whole-cohort decompression.

Target artifacts contain ordered window IDs, chromosomes, starts/ends, positions, IN/OUT class, legacy annotation IDs, GC, mappability, reference bases, chromosome/arm metadata, and source checksums. Preserve original coordinates as explicit fields; any zero-based internal view needs a tested reversible conversion.

Prepared artifacts contain per-sample integer counts, normalized values, diagnostic summaries, target checksum, and effective filter settings. Results contain HSLM tables, FastCall tables, both VCF forms, plots, and a run manifest. Preserve the original result basenames, including `HSLMResults_<sample>.txt`, `FastCallResults_<sample>.txt`, `EXCAVATORWindowCall_<sample>.vcf`, and `EXCAVATORRegionCall_<sample>.vcf`, beneath `Results/<sample>/`.

Manifests record schema version, compatibility profile, source commit/package version, complete parameters, input hashes, reference identity, CPU dispatch path, compiler/build identity, thread count, and completion status. Reject mixing counts from a different target, even if array lengths happen to match.

Provide an optional development converter from legacy RData to the canonical typed schema. Compare imported objects, not RData bytes. R is permitted for this converter/oracle only; it must not become a hidden production dependency.

## 5. Port each computational stage

### 5.1 Target generation and reference features

1. Export golden target rows from `FilterTarget.R`, retaining exact order and generated IDs.
2. Port chromosome selection, flank handling, off-target window construction, sorting, and gap filtering into readable Python/NumPy functions, preserving the literal legacy rules.
3. Test start/end equality, zero eligible gaps, chromosomes without targets, overlapping and unsorted target rows, naming prefixes, and windows spanning gap boundaries.
4. Implement Python feature extraction using existing FASTA/BigWig readers, preserving ambiguous-base handling, missing coverage, denominators, selected track fields, and reference-base offsets. Verify exact rather than approximate track summaries against the oracle.
5. Compare every target row and reference base exactly; compare GC/mappability with documented precision limits.
6. Only then replace repeated scans with interval indexing, sweeps, or cached sequence blocks, while keeping the legacy predicates intact.

Acceptance: the new target artifacts decode to the same ordered data used by preparation and VCF generation. A target mismatch blocks downstream parity claims.

### 5.2 Read filtering and counting

1. Capture the legacy selected `(chromosome, position)` stream on tiny BAM fixtures and the acceptance sample.
2. Test each relevant flag independently, MAPQ 0 and threshold-adjacent values, supplementary/secondary alignments, paired mates, duplicate positions, and contig boundaries. Reproduce what the code retains; do not add conventional read-quality policies.
3. Translate the `MakeReadCount.R`/`EXOMECOUNT` state machine together into a clear Python reference, including chunk boundaries and residual/final interval behavior. Use equivalent batched NumPy operations where proven safe.
4. Use indexed or sequential traversal through the Python reader, selected by measurement. Indexed and sequential paths must return identical counts in the required order. CRAM requires the correct explicit reference and independent parity fixtures. Explicitly adapt reader coordinates to the legacy SAM-position convention.
5. Keep contiguous count arrays and eliminate temporary text streams after equivalence is established. Profile decoding, Python record extraction, and counting separately. Add a small batch C++ counting kernel only if counting remains costly; it must accept position arrays and explicit residual state, leaving I/O and policy in Python. If extraction itself dominates, evaluate a narrowly scoped native-backed batch reader before expanding custom C++.
6. Test 499,999 / 500,000 / 500,001 selected-read boundaries, reads before the first interval, exactly on endpoints, gaps, final-window-only reads, and empty streams.

Acceptance: exact integer counts at every window. A one-read discrepancy is a failure, even if the final calls currently agree.

### 5.3 Normalization and experimental design

Implement reference-compatible median, quantile, variance, missing-value, and bin-boundary helpers in Python/NumPy first. Keep normalization and experimental-design arithmetic in Python modules; use a small compatibility helper only if a demonstrated mismatch requires it. Port size → MAP → GC correction in precisely that order, with separate IN/OUT behavior and zero replacement. Record every intermediate vector.

Reproduce the R coercion path as well as the arithmetic: `cbind` with character metadata creates character matrices, followed later by numeric parsing. Export values immediately before and after these boundaries to detect effective rounding. A port retaining more precision can still change a call.

Pool normalized controls in the original order and then divide by their count. Compute sample/reference ratios and separate IN/OUT median centering exactly as the legacy inference path. Parallelize independent samples only after each sample's arithmetic is stable.

Acceptance: exact masks and row order, numerical agreement at every correction stage, and no changes to later segmentation/calls. No LOWESS, pseudocount, robust estimator, or revised pooling rule is introduced here.

### 5.4 HSLM scalar reference port

Translate the R wrapper into Python and the demanding Fortran numerical core into C++. Python owns parameter estimation, state-grid construction, distance covariates, arm splitting, and segment reconstruction/filtering. C++ owns transition/emission evaluation and Viterbi/traceback. Keep a small-input Python reference for the latter so the complete algorithm can be understood and stepped through without reading C++.

* Estimate `mi`, `smu`, and `sepsilon` from the same globally selected data and quantiles. Match R's sample-variance convention and quantile interpolation.
* Reproduce the single-profile state grid including binary64 values, not just nominal decimal labels. Multi-profile helper support can be deferred if it is not reachable from the supported CLI; record that scope explicitly.
* Build distance-dependent transition probabilities using the existing expression:

  `eta_i = theta + (1 - theta) * exp(log(theta) / (delta_position_i / D_norm))`.

* Match Fortran's transition orientation, column-major indexing, emission formula, PI constant, and `ELNSUM` evaluation order. Substituting `log1p` may alter rounding and belongs after baseline proof.
* Translate Viterbi with binary64 scores and the original strict `>` comparison, which retains the first state on exact finite ties. Test the terminal selection and traceback independently.
* Preserve position integer coercion, centromere split rules, omitted/failed arm handling, and global parameter scope.
* Reproduce `SortState`, `FilterSeg`, and median segment reconstruction, including single-window and short-segment behavior.

The Fortran Viterbi routine declares a smaller transition shape than the flattened storage it accesses. Model the intended observed storage safely in C++; do not translate out-of-bounds declarations literally. Check small exported matrices and paths against the working reference build.

First deliver a transparent O(T K²) implementation, with K = 21 in the present single-profile CLI. Expose test-only traces for transition blocks, emissions, dynamic-programming scores, predecessors, paths, breakpoints, and filtered segment values.

After parity, reduce memory: generate each transition block on demand, retain two score vectors, and store traceback in a suitably sized integer buffer. This replaces O(T K²) transition storage with O(K²), while traceback remains O(T K). Compare against the transparent implementation on every fixture.

There is a possible later O(T K) recurrence because off-diagonal transitions share destination-dependent values. Treat this as a separate experiment: selecting top predecessors and reassociating sums can alter ties or rounding. It is not a prerequisite for the first release and must pass the same decision gate.

### 5.5 FastCall scalar reference port

Keep the actual algorithm visible in `fastcall.py`, with a readable Python/NumPy reference for accelerated calculations. Python owns steps 1–3, the iteration loop and stopping policy in step 6, assignment policy in step 7, and all of step 8. C++ batch kernels perform steps 4–5 and posterior evaluation in steps 6–7, using bounds and parameters supplied by Python. Do not replace this with a standard Gaussian-mixture estimator:

1. Construct the same segment summary from consecutive segment values, including any boundary quirks. Feed one value per resulting segment, without adding window-count weighting.
2. Apply cellularity correction to the same input with the same `2^-5` floor when cellularity is below one. Preserve the distinction between corrected calling values and uncorrected summary values used for reported CN.
3. Initialize fixed means `[-3, -1, 0, 0.58, 1]`, deviations, range-dependent deviation estimates, and priors exactly as `StartCond`/`EMFastCall` do.
4. Implement truncated Gaussian E-step terms using compatible normal density/CDF calculations. Preserve inclusive bounds, zero-denominator nearest-component fallback, first-nearest ties, and the explicit positive-infinity replacement.
5. Update deviations and priors with fixed means, including empty-component prior fallback and tiny-deviation handling. Preserve reduction order and prior renormalization.
6. Reproduce the stopping statistic, R matrix/vector recycling in its multiplication, the `1e-5` threshold, and the 1,000-iteration cap. A conventional log-likelihood convergence test is a later algorithm change.
7. Compute final untruncated posteriors as the original does, then assign `-2, -1, 0, 1, 2` labels using characterized tie behavior.
8. Reproduce significant-segment filtering, fractional CN, integer rounding, and output probabilities in the text and VCF writers separately.

For difficult CDF/tail behavior, compare a candidate native implementation against R across central and extreme inputs before selecting it. If reference numerical routines are reused, retain their notices and review redistribution requirements. Merely agreeing with another modern statistics library is insufficient.

Trace initial parameters, each iteration's deviations/priors/stopping statistic, final posterior matrix, and assignments. Tests must cover every truncation boundary, empty components, all-normal input, one segment, extreme values, underflow, exact/near ties, and cellularity below one.

Acceptance: same iteration behavior on deterministic characterization cases and exactly the same segment labels and reported discrete values on all acceptance fixtures. Continuous values obey the numerical contract.

## 6. SIMD, AVX, memory, and parallel execution

SIMD means “single instruction, multiple data”; AVX is an x86 SIMD instruction family. Treat SIMD as an optional implementation strategy, not a required machine capability.

First establish readable Python reference kernels against the legacy oracle, then qualify scalar C++ kernels against both. The native compatibility backend uses strict floating-point settings. Disable fast-math and uncontrolled FP contraction in compatibility builds. Inspect generated code/compiler reports where necessary: source-level scalar loops alone do not guarantee a particular floating-point execution order.

| Kernel | Candidate optimization | Main compatibility constraint |
|---|---|---|
| GC/base feature extraction | Python reader block access and NumPy operations; no planned custom SIMD | ambiguous-base denominator and endpoints |
| Filtering/counting | existing native-backed readers, batching; optional measured counting kernel | exact selected-read stream and legacy residual behavior |
| Normalization | Python/NumPy bin memberships, medians, vector transforms | bin equality rules, median definition, missing values |
| HSLM emissions | vector arithmetic across windows/states | identical constants and expression order |
| HSLM predecessor search | vector candidates with stable index selection | strict tie policy and floating-point association |
| FastCall density/E-step | structure-of-arrays, SIMD across segments | exp/CDF accuracy, underflow, bounds |
| FastCall M-step | vector element products; ordered accumulation | reduction order and stopping threshold |
| Ratios and CN transforms | Python/NumPy elementwise operations | log/exp and output rounding |

Prioritize elimination of subprocesses, text intermediates, repeated allocation, and huge transition tensors before hand-written intrinsics. Profile real stage time first: BAM decompression or target I/O can dominate total runtime even when segmentation accelerates substantially.

Compile portable scalar code plus optional AVX2 variants in isolated translation units. Add AVX-512 only if representative hardware measurements justify it. Use runtime CPU **and OS-state** detection before dispatch, and keep a forced-scalar mode for diagnosis. Do not build distributable wheels with global `-march=native`. ARM64 must retain a portable path; NEON optimization is a later independent backend.

Compiler function multiversioning is one possible dispatch mechanism, but support must be handled per toolchain; a small explicit dispatcher can provide a common interface. See [GCC's official multiversioning documentation](https://gcc.gnu.org/onlinedocs/gcc/Function-Multiversioning.html).

For each optimized kernel, compare scalar and optimized outputs on all fixtures, threshold-adjacent inputs, and randomized finite cases. Approximate vector exp/log implementations are not enabled merely because their average error is small. If an optimized backend changes calls, keep it disabled for the compatibility profile.

Parallelize samples first. Parallel chromosome/arm jobs may follow global HSLM parameter estimation; assemble them in the original chromosome order and run sample-level FastCall on the same complete segment set. Do not run independently fitted FastCall models per chromosome.

Use one bounded scheduler. Avoid nesting sample workers, chromosome workers, HTSlib decoding threads, and numerical thread pools at full size. Preserve deterministic reduction and control-pooling order. Do not change RNG consumption by scheduling tied assignments concurrently. Test `--threads 1`, `2`, and `4` for identical decisions.

Record wall time, CPU time, peak RSS, reads/windows/segments processed, throughput, hardware, compiler, dispatch path, storage, and thread count. Benchmark cold and warm I/O separately, with plots timed separately. Set speed/memory budgets after the baseline is measured; do not promise an unsupported “10×” speedup.

## 7. Differential testing against the original

### 7.1 Freeze a real oracle

The existing CI executes the three stages but does not assert equality with golden results. It also refers to a mutable container tag and downloads large reference assets. Build the following before porting numerical code:

1. Pin the baseline commit, container digest, platform, OS packages, shell/awk, R packages, Fortran compiler/flags, samtools/bedtools/UCSC tools, locale, and relevant environment settings. Capture a resolved dependency manifest; a mutable tag is insufficient.
2. Record checksums and sizes for every BAM/BAI or CRAM/index, BED, FASTA/index, BigWig, chromosome/gap/centromere table, YAML file, and parameters file.
3. Run the original on the `.test` paired sample and, when identified, the user's exact acceptance input. The FASTA and BigWig referenced by `.test/config.yaml` are not currently in the checked-in fixture listing; fetching and hashing them is part of this milestone.
4. Repeat runs at fixed thread count, then at different counts. Determine whether outputs are deterministic before choosing tolerances or golden calls.
5. Preserve the untouched final artifacts. Export RData objects through a development-only R script into typed binary arrays and metadata without decimal precision loss.
6. Capture checkpoints using a separate instrumented oracle copy and verify that instrumentation does not alter final calls relative to the untouched run.

If the exact original fails on an input, record its failure separately. Do not patch it silently and call that the original baseline. If undefined behavior prevents a stable reference result, resolve and document the baseline limitation before claiming equivalence.

### 7.2 Handle nondeterminism explicitly

Capture R RNG kind/version and seed/state around assignments if ties affect calls. The preferred reproducibility harness sets a fixed seed in a documented oracle wrapper while leaving the algorithm unchanged; retain an unmodified-run comparison too. Reproducing an existing output may require its RNG state.

If no ties occur on the acceptance sample, prove that through recorded posterior margins. If tied calls vary and the historical seed is unavailable, report that exact run reconstruction is unresolved. Do not accept arbitrary matching sets of labels or silently switch to a deterministic tie policy. A deterministic policy is a later behavior change unless proven equivalent for the supported cases.

### 7.3 Checkpoint hierarchy

| Checkpoint | Exact comparisons | Numerical comparisons |
|---|---|---|
| Target | row IDs/order, coordinates, class, reference bases | GC, mappability |
| Alignment selection | selected positions/flags in small diagnostic fixtures | none |
| Counting | all per-window counts, masks | none |
| Normalization | bin memberships, row order, NA masks | each correction, normalized counts |
| Experimental design | sample/control mapping and order | pooled control, ratios, centered ratios |
| HSLM parameters | state count/order, arm membership | state values, variance, model parameters |
| HSLM core | paths, predecessors on diagnostic fixtures, breaks | transition/emission/scores |
| Segment filtering | retained breakpoints, membership | segment medians |
| FastCall | segment membership, labels, diagnostic iteration counts | priors, deviations, stopping statistic, posteriors |
| Final outputs | record order, intervals, labels, CN, REF/ALT, meaningful VCF fields | segment values, CNF, probabilities |

Compare arrays with an explicit field contract, for example `abs(a-b) <= atol + rtol*abs(reference)`, plus separate equality checks for NA/NaN locations and infinity signs. Measure observed drift before setting each tolerance; provisional exploration thresholds are not release thresholds. Never widen a tolerance merely to make a failing port pass. ULP reports and decision margins help diagnose discrepancies.

Canonicalize only known volatile metadata such as file date and run-directory provenance. Do not sort away incorrect record order, merge nearby intervals, use overlap percentages in place of exact breakpoints, or round probabilities until differences disappear. Test output rendering/rounding separately from full-precision arrays.

### 7.4 Test layers and fixtures

**Python unit tests:** configuration, interval predicates, counting reference/state transitions, normalization/statistics conventions, segment filtering, EM control/stopping, assignment policy, and rounding. **Native unit tests:** transition/emission blocks, Viterbi tie/traceback, Gaussian terms, E/M updates, and any justified counting kernel. Use tiny manually inspectable cases plus oracle exports.

**Differential tests:** compare legacy R/Fortran against the Python implementation, and compare each custom C++ kernel against both the Python reference and legacy checkpoints. Compare the first differing checkpoint, not only the final VCF. Include conversion boundaries and exceptional-value masks.

**Hybrid integration tests:** old target → new prepare; old prepared arrays → new analyze; new target → old prepare via development conversion; new prepared arrays → old analyze. Temporary converters localize errors without making R a production dependency.

**End-to-end acceptance:** execute all three new commands from BAM/CRAM and compare final HSLM/FastCall/VCF outputs to the pinned oracle. Exact interval/label agreement is mandatory, including the no-call/header-only case.

**Coverage matrix:** paired and pooling designs; one and multiple controls; different sample ordering; nondefault HSLM/FastCall parameters; cellularity below one; IN/OUT mixtures; chromosome prefixes; both centromere arms; short/empty inputs where the oracle succeeds; repeated positions; zero counts; all five calling states; threshold-adjacent values; scalar and each enabled SIMD path; thread counts 1/2/4.

Add synthetic fixtures specifically to exercise positive calls and all states. A shallow dataset producing no calls alone cannot validate FastCall. Add a larger real sequencing fixture for performance and realistic normalization distributions; access and storage arrangements may be external to Git.

**Package tests:** install the built wheel in a clean environment; execute the documented commands; verify no R/Perl/Fortran runtime is present or invoked. Test supported Linux x86-64 and macOS ARM64 builds independently, with cross-platform exact decisions and qualified numerical tolerances.

**Safety and robustness tests:** C++ address/undefined-behavior sanitizers, invalid shapes/strides at bindings, interrupted writes, concurrent output-directory access, spaces in paths, mismatched target manifests, and corrupted input handling. Reproducing legacy scientific behavior does not require reproducing unsafe memory access.

Use metamorphic checks only after proving the transformation is an invariant of the reference; changing read chunking or interval ordering may expose a legacy quirk rather than an error in the port.

### 7.5 Regression reporting and CI

Each mismatch report should contain sample, chromosome/arm, first differing window/segment, reference/new values, absolute/relative/ULP error, relevant decision margin, effective parameters, and a runnable minimal reproduction.

Run small native, binding, CLI, and differential fixtures on every PR. Run the mandatory acceptance dataset for every compatibility-affecting change and release. Run the large BAM/reference suite and SIMD hardware matrix on scheduled or dedicated workers. Cache assets by content hashes; verify downloaded content before use.

Keep goldens immutable by default. A PR changing expected calls must identify a deliberate post-parity behavior change; a general “regenerate snapshots” command must not conceal regressions.

The parity report must list baseline identity, fixture hashes, exact mismatch counts, numerical maxima, backend/thread/platform results, excluded volatile fields, and unresolved limitations. “The pipeline completed” is not acceptance evidence.

## 8. Implementation milestones and exit gates

Current implementation status: **M2–M4 are implemented for the locally qualified
BAM/single-profile scope**. Python preparation and paired/pooling analysis use a
converted original target. Counts and prepared matrices match the supplied BAMs;
new prepare → analyze preserves original calls. See [M4 evidence and limits](preparation.md).
Live CI now also runs original versus current preparation on small BAMs, in addition
to analysis-only comparisons. Cross-platform M4 package tests and live concordance passed in
[CI run 36043110098](https://github.com/danilotat/excavator2-py/actions/runs/36043110098).
**M5.3 is complete for captured fixtures**: Python geometry and reference features
match original outputs, including captured failure cases. pyBigWig and pysam
provide native-backed reference access. Live feature comparison is added to CI.
Full target generation is not yet integrated; M5.4 is CLI/artifact integration and
end-to-end acceptance. See [M5 contract and remaining steps](target-porting.md).

Estimates below are rough focused engineering effort for one experienced developer, not calendar commitments. They exclude reference download delays and expand if oracle instability or numerical incompatibility is substantial. Complete exit gates in dependency order.

| Milestone | Indicative effort | Deliverables and required exit gate |
|---|---|---|
| M0: behavior inventory and oracle | 1–2 weeks | Pinned legacy environment, hashed inputs, repeatability report, checkpoint exporter, quirk register; usable golden outputs for acceptance data |
| M1: package skeleton and artifact contract | 1 week | Buildable wheel, C++ extension, CLI argument/config tests, manifests, typed arrays, legacy converter; clean install and artifact roundtrip succeed |
| M2: Python FastCall and scalar kernels | 1–2 weeks | Readable Python algorithm/reference, batch C++ E/M/posterior kernels, and full FastCall trace; differential tests pass on exported segments, all states, boundaries and ties |
| M3: Python HSLM control and scalar core | 2–3 weeks | Python setup/filtering/reference, C++ numerical recurrence, Python paired/pooling ratios and writers; new analyze matches old prepared inputs exactly in decisions |
| M4: preparation | 2–3 weeks | Python readers/filter policy, batched counts, NumPy normalization; optional counting kernel only after profiling; exact counts, matching normalization checkpoints, new prepare → analyze matches oracle |
| M5: target generation | 1–2 weeks | Python geometry/reference features and artifact resolver; all-new target → prepare → analyze matches acceptance calls |
| M6: full compatibility qualification | 1–2 weeks | Complete differential matrix, documented quirks, plot/data checks, independent clean wheel run; signed-off parity report with zero unexplained decision mismatches |
| M7: performance work | 2–3 weeks | Measured memory improvements, bounded concurrency, optional SIMD; every enabled path passes existing parity gates and measured improvements justify complexity |
| M8: release qualification | 1 week | Reproducible source/wheel builds, notices, user migration guide, supported-platform matrix, published benchmark/parity reports |

The order intentionally delivers analysis from exported prepared data early, then replaces preparation and target generation. This limits the number of changing upstream variables while debugging HSLM and FastCall. Compatibility tests remain mandatory at every milestone, not a final testing phase.

A practical task breakdown within each milestone is: capture reference inputs/outputs; write the behavioral contract; implement readable Python/NumPy logic; run differential tests; add a small batch binding only for the designated or measured hot kernel; integrate the CLI stage; document discrepancies; then accept the milestone. Keep algorithm porting and optimization in separate reviewable changes.

## 9. Risks and decisions to resolve during implementation

| Risk | Resolution strategy |
|---|---|
| Original run is unstable or fails | Pin environment, reproduce, isolate seed/undefined behavior, document unsupported cases; do not invent golden outputs |
| Small floating-point drift changes breakpoints | Locate first differing score, match conversion/math/reduction behavior, retain strict scalar fallback |
| Generic libraries implement different statistics | Differential-test each helper; retain explicitly compatible implementations |
| Legacy counting bugs depend on chunk boundaries | Characterize the wrapper and state machine together; preserve observable successful behavior |
| A faster HSLM recurrence changes tie behavior | Keep O(T K²) reference and require independent proof/tests before enabling an alternative |
| New schemas change coordinate meaning | Preserve source coordinates and test every I/O adapter and reference-base lookup |
| Full human reference data is too large for PR CI | Hash-addressed external assets plus tiny exported fixtures; mandatory full release gate |
| SIMD helps kernels but not wall time | Profile complete stages and report I/O/plot costs separately |
| Dependency wheels are difficult to distribute | Start clean-wheel tests early; keep native dependency surface small and versions recorded |

The repository's `LICENSE` identifies Attribution-NonCommercial-ShareAlike 3.0 Unported. Preserve existing attribution and notices during the port; do not assume that rewriting code changes the project's license. Include dependency/license review as a release packaging task rather than choosing a new license in this roadmap.

Decisions that need evidence during M0/M1: the exact user acceptance dataset if different from `.test`; supported CPython/platform versions; pinned oracle toolchain; R tie/RNG behavior; per-field tolerances; Python-accessible BigWig adapter; legacy artifact import scope; and benchmark hardware/data. None requires changing the requested CLI or relaxing call equivalence.

## 10. After parity: controlled bug fixes

Freeze and tag the compatibility release and its oracle/fixture manifest first. Then address known bugs one at a time:

1. Add a minimal reproducer showing the legacy behavior.
2. Write the intended corrected behavior and scientific rationale separately from the compatibility test.
3. Implement the correction under a versioned behavior profile or clearly documented release change. Keep the compatibility profile available for reproducibility.
4. Compare both profiles across the complete corpus and report changed intervals, states, probabilities, and affected samples.
5. Revalidate performance and documentation without replacing the old goldens.

Likely investigation candidates include MAPQ filtering, terminal read counting, coordinate/reference-base handling, target gap filtering, cross-boundary segment grouping, missing/short arms, and FastCall tie/convergence behavior. These remain deferred investigations until the all-new pipeline produces the same calls on the acceptance input.

## 11. Definition of done for the first port

- The three requested commands work from a clean installed package with compatible YAML inputs.
- HSLM and FastCall hot kernels run in the small compiled core; other stages use readable Python and efficient native-backed libraries. Any additional custom C++ is justified by profiling or a demonstrated compatibility need. Production execution requires no legacy language runtime.
- A Python developer can follow parameter setup, transformations, filtering, convergence, and calling policy in Python. Every custom kernel has a small-input Python reference and documented inputs/outputs.
- Target identities, read counts, HSLM boundaries, call intervals, labels, and reported discrete values match the pinned reference on all mandatory datasets.
- Continuous values meet documented tolerances, and every enabled platform/backend/thread configuration passes the decision gate.
- The exact user-provided acceptance file is included once identified; matching only a substitute fixture is insufficient.
- Artifacts are versioned, traceable, restart-safe, and reject incompatible target/preparation combinations.
- Supported wheels work on machines without AVX; optimized kernels use verified dispatch and retain a scalar fallback.
- Compatibility quirks are documented, suspected bug fixes are deferred, and golden results have not been silently rewritten.
- A reproducible parity report and measured performance report accompany the release.
