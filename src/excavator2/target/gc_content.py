"""
GC content calculation for target windows.

This module provides functions for calculating GC content of analysis
windows using the reference genome FASTA file.
"""

from pathlib import Path
from typing import List, Union, Optional
import logging

import numpy as np

from excavator2.io.fasta import FastaReader
from excavator2.target.filter import AnalysisWindow

logger = logging.getLogger(__name__)


def calculate_gc_content(windows: List[AnalysisWindow], fasta_path: Union[str, Path]) -> np.ndarray:
    """Calculate GC content for analysis windows.

    Args:
        windows: List of analysis windows
        fasta_path: Path to reference genome FASTA file (must be indexed)

    Returns:
        Array of GC fractions (0.0-1.0), same order as input windows
    """
    if not windows:
        return np.array([], dtype=np.float64)

    logger.info(f"Calculating GC content for {len(windows)} windows")

    # Extract coordinates
    chromosomes = [w.chrom for w in windows]
    starts = [w.start for w in windows]
    ends = [w.end for w in windows]

    # Calculate GC content
    with FastaReader(fasta_path) as reader:
        result = reader.get_gc_content(chromosomes, starts, ends)

    logger.info(f"Mean GC content: {result.mean_gc:.2%}")
    return result.gc_content


def calculate_gc_content_by_chromosome(
    windows: List[AnalysisWindow], fasta_path: Union[str, Path]
) -> dict:
    """Calculate GC content grouped by chromosome.

    Args:
        windows: List of analysis windows
        fasta_path: Path to reference genome FASTA file

    Returns:
        Dictionary mapping chromosome to array of GC fractions
    """
    if not windows:
        return {}

    # Group windows by chromosome
    windows_by_chrom: dict = {}
    indices_by_chrom: dict = {}

    for i, window in enumerate(windows):
        if window.chrom not in windows_by_chrom:
            windows_by_chrom[window.chrom] = []
            indices_by_chrom[window.chrom] = []
        windows_by_chrom[window.chrom].append(window)
        indices_by_chrom[window.chrom].append(i)

    # Calculate GC content for each chromosome
    result = {}
    with FastaReader(fasta_path) as reader:
        for chrom in sorted(windows_by_chrom.keys()):
            chrom_windows = windows_by_chrom[chrom]

            chromosomes = [w.chrom for w in chrom_windows]
            starts = [w.start for w in chrom_windows]
            ends = [w.end for w in chrom_windows]

            gc_result = reader.get_gc_content(chromosomes, starts, ends)
            result[chrom] = gc_result.gc_content

            logger.info(
                f"  {chrom}: {len(chrom_windows)} windows, mean GC={np.mean(gc_result.gc_content):.2%}"
            )

    return result


def get_window_gc_stats(gc_content: np.ndarray, windows: List[AnalysisWindow]) -> dict:
    """Calculate GC content statistics.

    Args:
        gc_content: Array of GC fractions
        windows: List of analysis windows

    Returns:
        Dictionary with GC statistics
    """
    if len(gc_content) == 0:
        return {
            "mean": 0.0,
            "median": 0.0,
            "std": 0.0,
            "min": 0.0,
            "max": 0.0,
            "in_target_mean": 0.0,
            "out_target_mean": 0.0,
        }

    # Overall stats
    stats = {
        "mean": float(np.mean(gc_content)),
        "median": float(np.median(gc_content)),
        "std": float(np.std(gc_content)),
        "min": float(np.min(gc_content)),
        "max": float(np.max(gc_content)),
    }

    # Stats by region class
    in_mask = np.array([w.region_class == "IN" for w in windows])
    out_mask = ~in_mask

    if np.any(in_mask):
        stats["in_target_mean"] = float(np.mean(gc_content[in_mask]))
    else:
        stats["in_target_mean"] = 0.0

    if np.any(out_mask):
        stats["out_target_mean"] = float(np.mean(gc_content[out_mask]))
    else:
        stats["out_target_mean"] = 0.0

    return stats
