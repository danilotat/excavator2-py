"""
Mappability extraction for target windows.

This module provides functions for extracting mappability scores
from BigWig files for analysis windows.
"""

from pathlib import Path
from typing import List, Union
import logging

import numpy as np

from excavator2.io.bigwig import BigWigReader
from excavator2.target.filter import AnalysisWindow

logger = logging.getLogger(__name__)


def extract_mappability(windows: List[AnalysisWindow], bigwig_path: Union[str, Path]) -> np.ndarray:
    """Extract mappability values for analysis windows.

    Args:
        windows: List of analysis windows
        bigwig_path: Path to mappability BigWig file

    Returns:
        Array of mappability values (0.0-1.0), same order as input windows
    """
    if not windows:
        return np.array([], dtype=np.float64)

    logger.info(f"Extracting mappability for {len(windows)} windows")

    # Extract coordinates
    chromosomes = [w.chrom for w in windows]
    starts = [w.start for w in windows]
    ends = [w.end for w in windows]

    # Get mappability values
    reader = BigWigReader(bigwig_path)
    result = reader.get_mappability(chromosomes, starts, ends, stat="mean")

    # Log warnings for invalid values
    n_invalid = result.n_regions - result.n_valid
    if n_invalid > 0:
        logger.warning(f"{n_invalid} windows have invalid mappability values")

    logger.info(f"Mean mappability: {result.mean_mappability:.3f}")
    return result.values


def extract_mappability_by_chromosome(
    windows: List[AnalysisWindow], bigwig_path: Union[str, Path]
) -> dict:
    """Extract mappability values grouped by chromosome.

    Args:
        windows: List of analysis windows
        bigwig_path: Path to mappability BigWig file

    Returns:
        Dictionary mapping chromosome to array of mappability values
    """
    if not windows:
        return {}

    # Group windows by chromosome
    windows_by_chrom: dict = {}

    for window in windows:
        if window.chrom not in windows_by_chrom:
            windows_by_chrom[window.chrom] = []
        windows_by_chrom[window.chrom].append(window)

    # Extract mappability for each chromosome
    result = {}
    reader = BigWigReader(bigwig_path)

    for chrom in sorted(windows_by_chrom.keys()):
        chrom_windows = windows_by_chrom[chrom]

        chromosomes = [w.chrom for w in chrom_windows]
        starts = [w.start for w in chrom_windows]
        ends = [w.end for w in chrom_windows]

        map_result = reader.get_mappability(chromosomes, starts, ends, stat="mean")
        result[chrom] = map_result.values

        mean_map = (
            float(np.mean(map_result.values[map_result.valid_mask]))
            if map_result.n_valid > 0
            else 0.0
        )
        logger.info(f"  {chrom}: {len(chrom_windows)} windows, mean MAP={mean_map:.3f}")

    return result


def get_mappability_stats(mappability: np.ndarray, windows: List[AnalysisWindow]) -> dict:
    """Calculate mappability statistics.

    Args:
        mappability: Array of mappability values
        windows: List of analysis windows

    Returns:
        Dictionary with mappability statistics
    """
    if len(mappability) == 0:
        return {
            "mean": 0.0,
            "median": 0.0,
            "std": 0.0,
            "min": 0.0,
            "max": 0.0,
            "in_target_mean": 0.0,
            "out_target_mean": 0.0,
            "low_mappability_count": 0,
            "low_mappability_fraction": 0.0,
        }

    # Overall stats
    stats = {
        "mean": float(np.mean(mappability)),
        "median": float(np.median(mappability)),
        "std": float(np.std(mappability)),
        "min": float(np.min(mappability)),
        "max": float(np.max(mappability)),
    }

    # Stats by region class
    in_mask = np.array([w.region_class == "IN" for w in windows])
    out_mask = ~in_mask

    if np.any(in_mask):
        stats["in_target_mean"] = float(np.mean(mappability[in_mask]))
    else:
        stats["in_target_mean"] = 0.0

    if np.any(out_mask):
        stats["out_target_mean"] = float(np.mean(mappability[out_mask]))
    else:
        stats["out_target_mean"] = 0.0

    # Count low mappability regions (< 0.5)
    low_map_threshold = 0.5
    low_map_count = int(np.sum(mappability < low_map_threshold))
    stats["low_mappability_count"] = low_map_count
    stats["low_mappability_fraction"] = low_map_count / len(mappability)

    return stats


def filter_low_mappability_windows(
    windows: List[AnalysisWindow], mappability: np.ndarray, min_mappability: float = 0.1
) -> tuple:
    """Filter windows with very low mappability.

    Args:
        windows: List of analysis windows
        mappability: Array of mappability values
        min_mappability: Minimum mappability threshold

    Returns:
        Tuple of (filtered_windows, filtered_mappability, n_filtered)
    """
    mask = mappability >= min_mappability
    n_filtered = int(np.sum(~mask))

    filtered_windows = [w for w, m in zip(windows, mask) if m]
    filtered_mappability = mappability[mask]

    if n_filtered > 0:
        logger.info(f"Filtered {n_filtered} windows with mappability < {min_mappability}")

    return filtered_windows, filtered_mappability, n_filtered
