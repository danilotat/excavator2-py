"""Time BAM decoding/counting and Python/NumPy correction separately."""

import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np

from excavator2.artifacts import load_manifest
from excavator2.normalization import normalize
from excavator2.reads import count_bam

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("bam")
parser.add_argument("target", type=Path)
args = parser.parse_args()
manifest = load_manifest(args.target, "target")
with np.load(args.target / manifest["preparation"]) as archive:
    target, gc, mappability = (archive[k] for k in ["target", "gc", "mappability"])
start = time.perf_counter()
counts = count_bam(args.bam, target, manifest["chromosomes"])
count_time = time.perf_counter() - start
start = time.perf_counter()
normalize(counts, target, gc, mappability)
normalization_time = time.perf_counter() - start
print(
    json.dumps(
        {
            "platform": platform.platform(),
            "windows": len(target),
            "bam_decode_and_count_seconds": count_time,
            "normalization_seconds": normalization_time,
            "threads": 1,
        },
        indent=2,
    )
)
