# M7 performance work

## M7.1: target-row selection

The first profile used the supplied 281,567-window acceptance dataset on an Apple
M4 Mac (16 GiB RAM), Python 3.12.6, NumPy 2.5.3, pysam 0.24.1 and pyBigWig 0.3.26.
Instrumented stage times were approximately 32.50 s for target, 5.12 s for
preparation and 3.61 s for analysis. These locate work; profiling overhead makes
them unsuitable for a speedup claim. The BAMs are shallow acceptance samples,
not a representative large cohort.

Target features consumed 28.19 s of the target profile. Exact BigWig means took
14.71 s; GC string counts took 4.94 s. Chromosome selection performed 6,476,042
whole-row regex searches. The M7.1 change indexes word-token chromosome occurrences
in a single pass, preserving original row order and matches in **any** column.
Chromosome names containing punctuation retain literal regex boundary matching.
Repeated occurrences still select a row once. This intentionally preserves the
original bare-chromosome selection problem; it is not a scientific correction.

No BigWig calculation, floating-point reduction, GC conversion, native kernel,
artifact schema or concurrency policy changes. This avoids introducing a new
numerical backend for a modest reduction in Python overhead.

## Measured result

Three fresh-process runs per implementation on the same local inputs, with no
cache eviction, include geometry, reference reads/hashes and artifact publication.
Every run compares all target arrays exactly with the qualified original target.

| Metric (median) | M6 baseline | Indexed selection |
| --- | ---: | ---: |
| Wall time | 30.340 s | 28.863 s |
| CPU time | 27.283 s | 25.759 s |
| Peak RSS | 475,676,672 bytes | 478,789,632 bytes |

Wall time is **4.9% lower (1.05× speedup)** in this measurement. Peak RSS is slightly
higher (about 3 MiB), so no memory improvement is claimed. There is no claim of
cold-disk, large-cohort or whole-pipeline speedup. The baseline used an isolated
M6 wheel and the candidate used the editable package with matching dependency
versions. The candidate wheel is qualified independently for correctness.

The reproducible runner is `benchmarks/target.py`; commands and caveats are in
`benchmarks/README.md`. Per-run wall/CPU/RSS measurements, input/source hashes and
runtime metadata are recorded in `benchmarks/reports/m7-target-index.json`.

## Remaining M7 work

BigWig queries remain the main measured target cost. Further batching must first
prove the same covered-base means and original decimal roundtrips on the M6
corpus. New native/SIMD code is not justified by the current HSLM/FastCall share.
Measure larger BAM/cohort workloads before changing preparation scheduling, and
measure memory independently before claiming a memory improvement. Existing
bounded sample preparation and single-worker analysis remain unchanged.

M7.1 is the first measured optimization; it does not complete the broader M7
performance milestone. Plotting remains deferred.

## M7.1 validation

All 189 tests pass in both the editable environment and an isolated wheel.
The optimized wheel's fresh supplied-data target → prepare → analyze pipeline
matches all 281,567 target windows and every prepared array against the qualified
baseline; scientific outputs match the saved original under the M6 gates.
The whole-row indexing regression includes bare names, cross-column matches,
repeated tokens, punctuation, underscores and Unicode word boundaries. Existing
BigWig threshold, precision and malformed-reference cases also pass. All seven
[remote CI jobs](https://github.com/danilotat/excavator2-py/actions/runs/36242892970)
passed for `2e75e61`, including Linux/macOS Python 3.11–3.13, installed wheels,
native sanitizers and the complete live legacy concordance chain. M7.1 is
complete; details are recorded in the benchmark report.
