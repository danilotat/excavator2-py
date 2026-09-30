"""Compare native and reference counting on a reproducible synthetic BAM."""

import json
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import median

import numpy as np
import pysam

from excavator2 import _core
from excavator2.reads import RegionCounter, count_bam
from excavator2.reference.bam_counts import count_positions, selected_chunks


def main():
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "reads.bam"
        n = 1_000_000
        with pysam.AlignmentFile(
            path, "wb", header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": n + 100}]}
        ) as bam:
            read = pysam.AlignedSegment()
            read.query_name = "benchmark"
            read.reference_id = 0
            read.mapping_quality = 60
            read.query_sequence = "A" * 50
            read.cigarstring = "50M"
            for i in range(n):
                read.reference_start = i
                bam.write(read)
        pysam.index(str(path))
        starts = np.arange(1, n, 50, dtype=np.int64)
        ends = starts + 100

        def reference():
            with pysam.AlignmentFile(path, "rb") as bam:
                return count_positions(selected_chunks(bam, "chr1"), starts, ends)

        plan = _core.RegionPlan("chr1", starts, ends)
        counter = RegionCounter()

        def native():
            return counter(path, (None, plan), 20)

        expected = reference()
        np.testing.assert_array_equal(native(), expected)
        timings = {}
        for name, operation in [("python", reference), ("cpp", native)]:
            samples = []
            for _ in range(5):
                start = time.perf_counter()
                operation()
                samples.append(time.perf_counter() - start)
            timings[name] = median(samples)
        for workers in [1, 4]:
            start = time.perf_counter()
            with ThreadPoolExecutor(max_workers=workers) as pool:
                for result in pool.map(lambda _: native(), range(8)):
                    np.testing.assert_array_equal(result, expected)
            timings[f"cpp_8_samples_{workers}_workers"] = time.perf_counter() - start
        target = np.column_stack([np.repeat("chr1", len(starts)), starts, ends])
        for workers in [1, 4, 16]:
            samples = []
            for _ in range(5):
                start = time.perf_counter()
                result = count_bam(path, target, ["chr1"], threads=workers)
                samples.append(time.perf_counter() - start)
                np.testing.assert_array_equal(result, expected)
            timings[f"single_bam_regions_{workers}_workers"] = median(samples)
        print(
            json.dumps(
                {
                    "reads": n,
                    "windows": len(starts),
                    "seconds": timings,
                    "speedup": timings["python"] / timings["cpp"],
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
