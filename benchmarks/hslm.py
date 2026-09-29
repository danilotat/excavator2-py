"""Measure single-arm HSLM including Python control; compare paths before timing."""

import json
import platform
import time

import numpy as np

from excavator2.hslm import estimate_parameters, segment

values = np.sin(np.arange(501) / 20)
positions = np.arange(501) * 10000
parameters = estimate_parameters(values)
results = {}
paths = {}
for backend in ["python", "native"]:
    times = []
    for _ in range(3):
        start = time.perf_counter()
        result = segment(
            values, positions, parameters, backend=backend, theta=0.0001, distance=100000
        )
        times.append(time.perf_counter() - start)
    paths[backend] = result.path
    results[backend] = float(np.median(times))
np.testing.assert_array_equal(paths["python"], paths["native"])
print(
    json.dumps(
        {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "windows": len(values),
            "states": 21,
            "repeats": 3,
            "median_seconds": results,
            "speedup": results["python"] / results["native"],
        },
        indent=2,
    )
)
