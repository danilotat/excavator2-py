# Benchmarks

After installing the package, run from the repository root:

```sh
python benchmarks/fastcall.py --segments 10003 --repeats 5
```

This deterministic synthetic benchmark checks posterior and iteration agreement
before timing both backends. JSON output records platform, Python, NumPy, input
size, iteration count, median full-fit times and speedup. It includes Python
control but excludes label assignment, HSLM and I/O; it does not measure pipeline
speedup or peak memory. Each fit runs scalar kernels on the calling thread.
Record compiler/build details separately when comparing environments.

Run `python benchmarks/hslm.py` for the scalar HSLM/reference comparison.
It checks paths and reports three-repeat medians for 501 windows and 21 states.

Run `python benchmarks/preparation.py BAM TARGET` to time BAM decoding/counting
and normalization separately against a converted target with preparation features.

## Complete target-stage timing and memory

```sh
/path/to/baseline/python benchmarks/target.py --config config.yaml \
  --expected-target qualified-target/ --output before/ --repeats 3
/path/to/candidate/python benchmarks/target.py --config config.yaml \
  --expected-target qualified-target/ --output after/ --repeats 3
```

Run the two commands sequentially on the same machine and inputs. Every repeat
uses a fresh process and output directory, records wall/CPU time and peak RSS,
and requires exact target-artifact equality against the qualified target. Timing
includes geometry, feature I/O, reference hashing and artifact publication; the
comparison itself runs after timing and RSS capture. JSON reports retain source
and input hashes, throughput, runtime and individual runs as well as medians.
Directories must not already exist. This script uses POSIX resource accounting
and supports the qualified Linux/macOS platforms.

No filesystem caches are evicted: report these as repeated/warm-cache local-file
measurements. Do not infer cold-storage performance, preparation/analysis speedup,
or large-cohort throughput from this target-stage result. For initial localization,
run a stage with `python -m cProfile -o stage.prof -m excavator2 ...`; profiling
overhead means those numbers are diagnostic rather than benchmark results.

## Original versus port in GitHub CI

The `Python port` workflow checks full scientific concordance before timing
comparable stages on its Ubuntu runner. The README table is generated from
`benchmark/report.json`; the latest published report is retained at
`benchmarks/reports/ci-latest.json`, and each run uploads its report and logs.

To reproduce locally with the installed port and Docker:

```sh
python tools/oracle/end_to_end.py --output .oracle/version-comparison
python tools/oracle/benchmark_versions.py --concordance .oracle/version-comparison --repeats 3
```

Use a fresh output directory. The original version uses the pinned source and
container from the concordance harness. Both versions receive the same synthetic
inputs: 4,113 windows, three prepared samples, and one paired analysis test sample.
Target timing uses one invocation per version from the full-pipeline comparison.
Paired analysis uses three measured invocations per version after an excluded
warm-up, alternating execution order and reporting medians. Scientific outputs
are checked after every analysis pair; differences fail the job.

Timing includes interpreter/process startup and artifact I/O. It excludes Docker
startup, image pulls, fixture generation, exports and comparison. Preparation is
not timed because the original preparation also generates plots, which the port
does not yet implement. Neither timed stage generates plots. Analysis uses one
worker and single-threaded native-library pools. Filesystem caches are not evicted.
The original runs inside Docker while the port runs on the host; on non-x86
machines Docker may additionally emulate the original's architecture. Local
numbers from such machines are not comparable to native x86 CI measurements.

These are small synthetic stage timings, not an end-to-end or real-WES speedup
claim. Hosted-runner noise is expected; performance ratios do not gate CI.
Correctness and valid timing records do. Compare the recorded revision, inputs
and environment before interpreting changes across runs.

Only successful pushes to this repository's `dev/porting` branch publish the
README table, after all package and concordance checks pass. Pull requests retain
read-only permissions. Publication validates the report's source revision and
updates only the marked README block plus the JSON report. A newer branch commit
causes publication to be skipped so an older result cannot overwrite newer work.
