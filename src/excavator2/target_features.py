"""Reference features matching TargetCreate.sh and its decimal text roundtrips."""

import re
from pathlib import Path

import numpy as np
import pyBigWig
import pysam


def target_features(target: np.ndarray, fasta: Path, bigwig: Path) -> dict:
    """Return per-chromosome gc, mappability and two-column first_base arrays.

    Uses indexed FASTA and exact covered-base BigWig means. Deliberately retains
    legacy grep selection and skipped FASTA rows, so feature lengths can differ.
    This extracts features only; it does not write target artifacts.
    """
    target = np.asarray(target, dtype=str)
    if target.ndim != 2 or target.shape[1] != 5 or not len(target):
        raise ValueError("target must be a nonempty five-column matrix")
    chromosomes = list(dict.fromkeys(target[:, 0]))
    prefix = "" if any("chr" in c for c in chromosomes) else "chr"
    lines = ["\t".join(row) for row in target]
    result = {}
    with pysam.FastaFile(str(fasta)) as reference, pyBigWig.open(str(bigwig)) as track:
        if not track.isBigWig():
            raise ValueError("mappability reference must be a BigWig file")
        lengths = dict(zip(reference.references, reference.lengths, strict=True))
        bw_lengths = track.chroms()
        for chromosome in chromosomes:
            # TargetCreate.sh uses grep -w on the entire row, not column one.
            pattern = re.compile(rf"(?<!\w){re.escape(chromosome)}(?!\w)")
            selected = [
                row for row, line in zip(target, lines, strict=True) if pattern.search(line)
            ]
            gc, mappability, first_base = [], [], []
            for chrom, start, end, _, _ in selected:
                chrom = prefix + chrom
                start, end = int(start), int(end)
                left = start - 1
                if left < 0 or end <= left:
                    raise ValueError(
                        "legacy feature extraction requires positive valid coordinates"
                    )
                # UCSC's sixth output column averages covered bases; no coverage is zero.
                bw_end = min(end, bw_lengths.get(chrom, 0))
                mean = None
                if left < bw_end:
                    mean = track.stats(chrom, left, bw_end, type="mean", exact=True)[0]
                mappability.append(float(format(0.0 if mean is None else mean, ".6g")))
                if chrom in lengths and end <= lengths[chrom]:
                    sequence = reference.fetch(chrom, left, end).upper()
                    fraction = (sequence.count("G") + sequence.count("C")) / len(sequence)
                    gc.append(float(format(fraction, ".6f")))
                # FRB is one base to the right of the GC/MAP left endpoint.
                if chrom in lengths and start < lengths[chrom]:
                    first_base.append((str(start), reference.fetch(chrom, start, start + 1)))
            if not first_base:
                raise ValueError(f"legacy first-base extraction has no rows for {chromosome}")
            result[chromosome] = {
                "gc": np.asarray(gc, dtype=float),
                "mappability": np.asarray(mappability, dtype=float),
                "first_base": np.asarray(first_base, dtype=str),
            }
    return result
