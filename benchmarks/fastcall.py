"""Measure both FastCall backends, with a correctness check before timing."""

import argparse
import json
import platform
import statistics
import time

import numpy as np

from excavator2.fastcall import fit_fastcall


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segments", type=int, default=10003)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if args.segments < 1 or args.repeats < 1:
        parser.error("segments and repeats must be positive")
    rng = np.random.default_rng(731)
    means = np.array([-3.0, -1.0, 0.0, 0.58, 1.0])
    values = rng.choice(means, args.segments) + rng.normal(0, 0.1, args.segments)
    expected = fit_fastcall(values, backend="python")
    actual = fit_fastcall(values, backend="native")
    np.testing.assert_allclose(actual.posterior, expected.posterior, rtol=1e-13, atol=2e-15)
    assert actual.iterations == expected.iterations
    timings = {}
    for backend in ("python", "native"):
        samples = []
        for _ in range(args.repeats):
            start = time.perf_counter()
            fit_fastcall(values, backend=backend)
            samples.append(time.perf_counter() - start)
        timings[backend] = statistics.median(samples)
    print(
        json.dumps(
            {
                "platform": platform.platform(),
                "python": platform.python_version(),
                "numpy": np.__version__,
                "segments": args.segments,
                "repeats": args.repeats,
                "iterations": actual.iterations,
                "median_seconds": timings,
                "speedup": timings["python"] / timings["native"],
                "scope": "synthetic segment fit including Python control; excludes HSLM/I/O/labels",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
