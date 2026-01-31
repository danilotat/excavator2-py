"""
Quality control plots for EXCAVATOR2 normalization.

This module generates bias plots showing the relationship between
read counts and various features (size, mappability, GC content)
before and after normalization.

These plots help assess the effectiveness of the normalization
pipeline at removing systematic biases.
"""

from pathlib import Path
from typing import Optional, Union
import logging

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from excavator2.prepare.readcount import SampleReadCounts
from excavator2.prepare.normalize import NormalizationResult
from excavator2.report.utils import (
    FIGURE_DPI,
    REGION_COLORS,
    compute_bin_statistics,
)

logger = logging.getLogger(__name__)


def create_qc_plots(
    raw_data: SampleReadCounts,
    norm_result: NormalizationResult,
    output_dir: Union[str, Path],
    format: str = "pdf",
) -> None:
    """Generate all QC plots for a sample.

    Creates bias plots showing raw vs normalized read counts
    against size, mappability, and GC content.

    Args:
        raw_data: Raw read count data
        norm_result: Normalized result
        output_dir: Directory to save plots
        format: Output format ('pdf', 'png', 'svg')
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_name = raw_data.sample_name
    logger.info(f"Creating QC plots for {sample_name}")

    # Extract data arrays
    raw_counts = raw_data.raw_counts
    norm_counts = norm_result.normalized_counts
    windows = raw_data.windows

    lengths = np.array([w.length for w in windows])
    mappability = np.array([w.mappability * 100 for w in windows])  # %
    gc_content = np.array([w.gc_content * 100 for w in windows])  # %
    in_mask = raw_data.get_in_target_mask()
    out_mask = ~in_mask

    # Generate plots
    # 1. Size bias (IN-target only)
    if np.any(in_mask):
        plot_size_bias(
            lengths[in_mask],
            raw_counts[in_mask],
            norm_counts[in_mask],
            output_path=output_dir / f"InSizeBias.{format}",
            title=f"{sample_name} - In-Target Size Bias",
        )
        logger.info(f"  Created InSizeBias.{format}")

    # 2. Mappability bias (IN-target)
    if np.any(in_mask):
        plot_mappability_bias(
            mappability[in_mask],
            raw_counts[in_mask],
            norm_counts[in_mask],
            output_path=output_dir / f"InMAPBias.{format}",
            title=f"{sample_name} - In-Target Mappability Bias",
        )
        logger.info(f"  Created InMAPBias.{format}")

    # 3. Mappability bias (OUT-target)
    if np.any(out_mask):
        plot_mappability_bias(
            mappability[out_mask],
            raw_counts[out_mask],
            norm_counts[out_mask],
            output_path=output_dir / f"OutMAPBias.{format}",
            title=f"{sample_name} - Off-Target Mappability Bias",
        )
        logger.info(f"  Created OutMAPBias.{format}")

    # 4. GC bias (IN-target)
    if np.any(in_mask):
        plot_gc_bias(
            gc_content[in_mask],
            raw_counts[in_mask],
            norm_counts[in_mask],
            output_path=output_dir / f"InGCBias.{format}",
            title=f"{sample_name} - In-Target GC Content Bias",
        )
        logger.info(f"  Created InGCBias.{format}")

    # 5. GC bias (OUT-target)
    if np.any(out_mask):
        plot_gc_bias(
            gc_content[out_mask],
            raw_counts[out_mask],
            norm_counts[out_mask],
            output_path=output_dir / f"OutGCBias.{format}",
            title=f"{sample_name} - Off-Target GC Content Bias",
        )
        logger.info(f"  Created OutGCBias.{format}")

    logger.info(f"QC plots saved to {output_dir}")


def plot_size_bias(
    lengths: np.ndarray,
    raw_counts: np.ndarray,
    norm_counts: np.ndarray,
    output_path: Optional[Union[str, Path]] = None,
    title: str = "Size Bias",
    bin_size: float = 5.0,
) -> plt.Figure:
    """Plot size/length bias showing raw vs normalized counts.

    Shows WMRC vs exon length with error bars indicating quartile ranges.

    Args:
        lengths: Window lengths in bp
        raw_counts: Raw read counts
        norm_counts: Normalized read counts
        output_path: Path to save figure (optional)
        title: Plot title
        bin_size: Bin size for length grouping (bp)

    Returns:
        Matplotlib Figure object
    """
    fig, axes = plt.subplots(2, 1, figsize=(15, 10), dpi=FIGURE_DPI)

    # Raw counts
    _plot_bias_panel(
        axes[0],
        lengths,
        raw_counts,
        bin_size,
        xlabel="Exon Length (bp)",
        ylabel="Raw WMRC",
        title=f"{title} - Raw",
        color="#1f77b4",
    )

    # Normalized counts
    _plot_bias_panel(
        axes[1],
        lengths,
        norm_counts,
        bin_size,
        xlabel="Exon Length (bp)",
        ylabel="Normalized WMRC",
        title=f"{title} - Normalized",
        color="#2ca02c",
    )

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)

    return fig


def plot_mappability_bias(
    mappability: np.ndarray,
    raw_counts: np.ndarray,
    norm_counts: np.ndarray,
    output_path: Optional[Union[str, Path]] = None,
    title: str = "Mappability Bias",
    bin_size: float = 5.0,
) -> plt.Figure:
    """Plot mappability bias showing raw vs normalized counts.

    Shows WMRC vs mappability percentage with error bars.

    Args:
        mappability: Mappability scores (0-100%)
        raw_counts: Raw read counts
        norm_counts: Normalized read counts
        output_path: Path to save figure (optional)
        title: Plot title
        bin_size: Bin size for mappability grouping (%)

    Returns:
        Matplotlib Figure object
    """
    fig, axes = plt.subplots(2, 1, figsize=(15, 10), dpi=FIGURE_DPI)

    # Raw counts
    _plot_bias_panel(
        axes[0],
        mappability,
        raw_counts,
        bin_size,
        xlabel="Mappability (%)",
        ylabel="Raw WMRC",
        title=f"{title} - Raw",
        color="#1f77b4",
    )

    # Normalized counts
    _plot_bias_panel(
        axes[1],
        mappability,
        norm_counts,
        bin_size,
        xlabel="Mappability (%)",
        ylabel="Normalized WMRC",
        title=f"{title} - Normalized",
        color="#2ca02c",
    )

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)

    return fig


def plot_gc_bias(
    gc_content: np.ndarray,
    raw_counts: np.ndarray,
    norm_counts: np.ndarray,
    output_path: Optional[Union[str, Path]] = None,
    title: str = "GC Content Bias",
    bin_size: float = 5.0,
) -> plt.Figure:
    """Plot GC content bias showing raw vs normalized counts.

    Shows WMRC vs GC percentage with error bars.

    Args:
        gc_content: GC content (0-100%)
        raw_counts: Raw read counts
        norm_counts: Normalized read counts
        output_path: Path to save figure (optional)
        title: Plot title
        bin_size: Bin size for GC grouping (%)

    Returns:
        Matplotlib Figure object
    """
    fig, axes = plt.subplots(2, 1, figsize=(15, 10), dpi=FIGURE_DPI)

    # Raw counts
    _plot_bias_panel(
        axes[0],
        gc_content,
        raw_counts,
        bin_size,
        xlabel="GC Content (%)",
        ylabel="Raw WMRC",
        title=f"{title} - Raw",
        color="#1f77b4",
    )

    # Normalized counts
    _plot_bias_panel(
        axes[1],
        gc_content,
        norm_counts,
        bin_size,
        xlabel="GC Content (%)",
        ylabel="Normalized WMRC",
        title=f"{title} - Normalized",
        color="#2ca02c",
    )

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)

    return fig


def _plot_bias_panel(
    ax: plt.Axes,
    feature: np.ndarray,
    counts: np.ndarray,
    bin_size: float,
    xlabel: str,
    ylabel: str,
    title: str,
    color: str = "#1f77b4",
) -> None:
    """Plot a single bias panel with error bars.

    Args:
        ax: Matplotlib axes
        feature: Feature values (x-axis)
        counts: Read counts (y-axis)
        bin_size: Bin size for grouping
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Panel title
        color: Plot color
    """
    # Compute binned statistics
    bin_centers, medians, q1, q3 = compute_bin_statistics(feature, counts, bin_size)

    if len(bin_centers) == 0:
        ax.text(
            0.5,
            0.5,
            "Insufficient data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14,
        )
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title(title, fontsize=14, fontweight="bold")
        return

    # Plot error bars (quartile range)
    errors = np.array([medians - q1, q3 - medians])
    ax.errorbar(
        bin_centers,
        medians,
        yerr=errors,
        fmt="o",
        color=color,
        ecolor=color,
        capsize=3,
        capthick=1.5,
        markersize=6,
        linewidth=1.5,
    )

    # Connect medians with line
    ax.plot(bin_centers, medians, "-", color=color, alpha=0.5, linewidth=1)

    # Add reference line at overall median
    overall_median = np.median(medians)
    ax.axhline(
        overall_median,
        color="grey",
        linestyle="--",
        alpha=0.7,
        linewidth=1,
        label=f"Median: {overall_median:.1f}",
    )

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(True, alpha=0.3)

    # Set reasonable y-limits
    ymax = np.max(q3) * 1.2
    ymin = max(0, np.min(q1) * 0.8)
    ax.set_ylim(ymin, ymax)


def create_combined_qc_plot(
    raw_data: SampleReadCounts,
    norm_result: NormalizationResult,
    output_path: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """Create a combined 2x3 QC plot showing all biases.

    This creates a compact summary of all bias corrections.

    Args:
        raw_data: Raw read count data
        norm_result: Normalized result
        output_path: Path to save figure (optional)

    Returns:
        Matplotlib Figure object
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), dpi=FIGURE_DPI)

    # Extract data
    raw_counts = raw_data.raw_counts
    norm_counts = norm_result.normalized_counts
    windows = raw_data.windows

    lengths = np.array([w.length for w in windows])
    mappability = np.array([w.mappability * 100 for w in windows])
    gc_content = np.array([w.gc_content * 100 for w in windows])
    in_mask = raw_data.get_in_target_mask()

    # Use IN-target data for this summary
    if np.any(in_mask):
        lengths = lengths[in_mask]
        mappability = mappability[in_mask]
        gc_content = gc_content[in_mask]
        raw_counts = raw_counts[in_mask]
        norm_counts = norm_counts[in_mask]

    # Row 1: Raw data
    _plot_bias_panel(
        axes[0, 0],
        lengths,
        raw_counts,
        5.0,
        "Exon Length (bp)",
        "Raw WMRC",
        "Size Bias - Raw",
        "#1f77b4",
    )
    _plot_bias_panel(
        axes[0, 1],
        mappability,
        raw_counts,
        5.0,
        "Mappability (%)",
        "Raw WMRC",
        "Mappability Bias - Raw",
        "#1f77b4",
    )
    _plot_bias_panel(
        axes[0, 2],
        gc_content,
        raw_counts,
        5.0,
        "GC Content (%)",
        "Raw WMRC",
        "GC Bias - Raw",
        "#1f77b4",
    )

    # Row 2: Normalized data
    _plot_bias_panel(
        axes[1, 0],
        lengths,
        norm_counts,
        5.0,
        "Exon Length (bp)",
        "Norm WMRC",
        "Size Bias - Normalized",
        "#2ca02c",
    )
    _plot_bias_panel(
        axes[1, 1],
        mappability,
        norm_counts,
        5.0,
        "Mappability (%)",
        "Norm WMRC",
        "Mappability Bias - Normalized",
        "#2ca02c",
    )
    _plot_bias_panel(
        axes[1, 2],
        gc_content,
        norm_counts,
        5.0,
        "GC Content (%)",
        "Norm WMRC",
        "GC Bias - Normalized",
        "#2ca02c",
    )

    fig.suptitle(f"QC Summary - {raw_data.sample_name}", fontsize=16, fontweight="bold")
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)

    return fig
