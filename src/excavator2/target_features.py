"""Reference features matching TargetCreate.sh and its decimal text roundtrips."""

import math
import re
from decimal import Decimal
from pathlib import Path

import numpy as np
import pyBigWig
import pysam


def legacy_decimal(text: str) -> float:
    """Replay pinned x86 R's decimal -> 64-bit significand -> binary64 rounding.

    R parses through extended precision. NumPy longdouble is only binary64 on
    macOS ARM, so use exact integer arithmetic for this small decimal boundary.
    This is for finite feature text, not a general replacement for R's parser.
    """
    numerator, denominator = Decimal(text).as_integer_ratio()
    if numerator == 0:
        return 0.0
    sign = -1 if numerator < 0 else 1
    numerator = abs(numerator)
    exponent = numerator.bit_length() - denominator.bit_length()
    if exponent >= 0:
        if numerator < denominator << exponent:
            exponent -= 1
    elif numerator << -exponent < denominator:
        exponent -= 1
    shift = 63 - exponent
    if shift >= 0:
        numerator <<= shift
    else:
        denominator <<= -shift
    significand, remainder = divmod(numerator, denominator)
    if 2 * remainder > denominator or (2 * remainder == denominator and significand % 2):
        significand += 1
    return sign * math.ldexp(float(significand), -shift)


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
                mappability.append(legacy_decimal(format(0.0 if mean is None else mean, ".6g")))
                if chrom in lengths and end <= lengths[chrom]:
                    sequence = reference.fetch(chrom, left, end).upper()
                    fraction = (sequence.count("G") + sequence.count("C")) / len(sequence)
                    gc.append(legacy_decimal(format(float(np.float32(fraction)), ".6f")))
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
