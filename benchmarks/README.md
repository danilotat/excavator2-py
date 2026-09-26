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
