"""
Log2 ratio computation for CNV analysis.

This module provides functions for computing log2 ratios between test and
control samples, including median centering and region-specific handling.
"""

from dataclasses import dataclass
from typing import List, Optional, Union
import logging

import numpy as np

from excavator2.prepare.normalize import NormalizationResult, load_normalized_counts
from excavator2.prepare.readcount import WindowData

logger = logging.getLogger(__name__)


@dataclass
class Log2RatioResult:
    """Result of log2 ratio computation.

    Attributes:
        sample_name: Test sample name
        control_name: Control sample name
        log2_ratios: Log2(test/control) values
        positions: Genomic positions for each window
        windows: Window metadata
        chromosomes: List of chromosomes
        in_target_mask: Boolean mask for IN-target windows
    """

    sample_name: str
    control_name: str
    log2_ratios: np.ndarray
    positions: np.ndarray
    windows: List[WindowData]
    chromosomes: List[str]
    in_target_mask: np.ndarray

    @property
    def n_windows(self) -> int:
        """Total number of windows."""
        return len(self.log2_ratios)

    def get_chromosome_data(self, chromosome: str) -> dict:
        """Get log2 ratio data for a specific chromosome.

        Args:
            chromosome: Chromosome name

        Returns:
            Dictionary with log2_ratios, positions, windows, in_target_mask
        """
        # Handle chr prefix variants
        chrom_variants = {chromosome}
        if chromosome.startswith("chr"):
            chrom_variants.add(chromosome[3:])
        else:
            chrom_variants.add(f"chr{chromosome}")

        mask = np.array([w.chrom in chrom_variants for w in self.windows])

        return {
            "log2_ratios": self.log2_ratios[mask],
            "positions": self.positions[mask],
            "windows": [w for w, m in zip(self.windows, mask) if m],
            "in_target_mask": self.in_target_mask[mask],
            "n_windows": int(np.sum(mask)),
        }

    def get_in_target_data(self) -> dict:
        """Get log2 ratio data for IN-target windows only.

        Returns:
            Dictionary with log2_ratios, positions, windows
        """
        mask = self.in_target_mask
        return {
            "log2_ratios": self.log2_ratios[mask],
            "positions": self.positions[mask],
            "windows": [w for w, m in zip(self.windows, mask) if m],
            "n_windows": int(np.sum(mask)),
        }

    def get_off_target_data(self) -> dict:
        """Get log2 ratio data for OFF-target windows only.

        Returns:
            Dictionary with log2_ratios, positions, windows
        """
        mask = ~self.in_target_mask
        return {
            "log2_ratios": self.log2_ratios[mask],
            "positions": self.positions[mask],
            "windows": [w for w, m in zip(self.windows, mask) if m],
            "n_windows": int(np.sum(mask)),
        }


def compute_log2_ratio(
    test_data: NormalizationResult,
    control_data: NormalizationResult,
    median_center: bool = True,
    separate_regions: bool = True,
    pseudocount: float = 1e-10,
) -> Log2RatioResult:
    """Compute log2 ratios between test and control samples.

    Calculates log2(test/control) for each window. Optionally performs
    median centering to remove global biases.

    Args:
        test_data: Normalized counts for test sample
        control_data: Normalized counts for control sample
        median_center: Center ratios by subtracting median (default: True)
        separate_regions: Center IN and OUT regions separately (default: True)
        pseudocount: Small value added to avoid log(0) (default: 1e-10)

    Returns:
        Log2RatioResult with computed ratios

    Raises:
        ValueError: If test and control have mismatched windows
    """
    logger.info(f"Computing log2 ratios: {test_data.sample_name} vs {control_data.sample_name}")

    # Validate window compatibility
    if len(test_data.windows) != len(control_data.windows):
        raise ValueError(
            f"Window count mismatch: test={len(test_data.windows)}, "
            f"control={len(control_data.windows)}"
        )

    # Validate chromosome order
    if test_data.chromosomes != control_data.chromosomes:
        logger.warning("Chromosome order differs between test and control")

    # Get counts
    test_counts = test_data.normalized_counts
    control_counts = control_data.normalized_counts

    # Add pseudocount to avoid division by zero
    test_safe = test_counts + pseudocount
    control_safe = control_counts + pseudocount

    # Compute log2 ratio
    log2_ratios = np.log2(test_safe / control_safe)

    # Build region mask
    in_target_mask = np.array([w.region_class == "IN" for w in test_data.windows])

    # Median centering
    if median_center:
        if separate_regions:
            # Center IN and OUT regions separately
            if np.any(in_target_mask):
                in_median = np.median(log2_ratios[in_target_mask])
                log2_ratios[in_target_mask] -= in_median
                logger.info(f"  IN-target median: {in_median:.4f} (centered)")

            out_mask = ~in_target_mask
            if np.any(out_mask):
                out_median = np.median(log2_ratios[out_mask])
                log2_ratios[out_mask] -= out_median
                logger.info(f"  OFF-target median: {out_median:.4f} (centered)")
        else:
            # Center all together
            global_median = np.median(log2_ratios)
            log2_ratios -= global_median
            logger.info(f"  Global median: {global_median:.4f} (centered)")

    # Extract positions
    positions = np.array([w.position for w in test_data.windows])

    return Log2RatioResult(
        sample_name=test_data.sample_name,
        control_name=control_data.sample_name,
        log2_ratios=log2_ratios,
        positions=positions,
        windows=test_data.windows,
        chromosomes=test_data.chromosomes,
        in_target_mask=in_target_mask,
    )


def compute_log2_ratio_pooled(
    test_data: NormalizationResult,
    control_samples: List[NormalizationResult],
    median_center: bool = True,
    separate_regions: bool = True,
    pseudocount: float = 1e-10,
) -> Log2RatioResult:
    """Compute log2 ratios against pooled control samples.

    For pooled experimental design, control samples are averaged
    to create a reference profile.

    Args:
        test_data: Normalized counts for test sample
        control_samples: List of normalized counts for control samples
        median_center: Center ratios by subtracting median (default: True)
        separate_regions: Center IN and OUT regions separately (default: True)
        pseudocount: Small value added to avoid log(0) (default: 1e-10)

    Returns:
        Log2RatioResult with computed ratios

    Raises:
        ValueError: If no control samples provided
    """
    if not control_samples:
        raise ValueError("At least one control sample required")

    logger.info(f"Creating pooled control from {len(control_samples)} samples")

    # Average control counts
    control_counts_list = [s.normalized_counts for s in control_samples]

    # Validate all have same length
    n_windows = len(test_data.windows)
    for i, counts in enumerate(control_counts_list):
        if len(counts) != n_windows:
            raise ValueError(f"Control sample {i} has {len(counts)} windows, expected {n_windows}")

    # Create pooled control (mean)
    pooled_counts = np.mean(np.array(control_counts_list), axis=0)

    # Create a synthetic NormalizationResult for the pooled control
    pooled_result = NormalizationResult(
        sample_name="PooledControl",
        normalized_counts=pooled_counts,
        windows=test_data.windows,  # Use test windows (should be identical)
        chromosomes=test_data.chromosomes,
        normalization_stats={},
    )

    return compute_log2_ratio(
        test_data,
        pooled_result,
        median_center=median_center,
        separate_regions=separate_regions,
        pseudocount=pseudocount,
    )


def apply_cellularity_correction(log2_ratios: np.ndarray, cellularity: float) -> np.ndarray:
    """Apply cellularity correction to log2 ratios.

    When tumor purity (cellularity) is less than 1.0, the observed
    log2 ratios are diluted by normal cell contamination. This function
    adjusts the ratios to estimate the true tumor copy number signal.

    Args:
        log2_ratios: Log2 ratio values
        cellularity: Tumor purity (0.0-1.0)

    Returns:
        Corrected log2 ratios
    """
    if cellularity >= 1.0:
        return log2_ratios.copy()

    if cellularity <= 0:
        raise ValueError(f"cellularity must be positive, got {cellularity}")

    # Convert log2 ratio to copy number fold change
    # R = 2^log2ratio
    # R_observed = cell * R_tumor + (1-cell) * 1
    # R_tumor = (R_observed - (1-cell)) / cell

    ratio = np.power(2, log2_ratios)
    corrected_ratio = (ratio - (1 - cellularity)) / cellularity

    # Clamp to avoid log of non-positive values
    min_ratio = 2 ** (-5)  # Minimum allowed ratio
    corrected_ratio = np.maximum(corrected_ratio, min_ratio)

    return np.log2(corrected_ratio)
