"""
Per-chromosome CNV visualization plots.

This module generates detailed plots for each chromosome showing:
1. Log2 ratio scatter plot with segmentation overlay
2. CNV calls as colored rectangles

These plots provide detailed visualization of the segmentation
and calling results for manual inspection and quality assessment.
"""

from pathlib import Path
from typing import List, Optional, Union
import logging

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

from excavator2.analyze.pipeline import AnalysisResult, CNVSegment, ChromosomeResult
from excavator2.analyze.ratio import Log2RatioResult
from excavator2.report.utils import (
    FIGURE_DPI,
    DEFAULT_FIGSIZE,
    SCATTER_SIZE,
    SCATTER_ALPHA,
    SEGMENT_COLOR,
    SEGMENT_LINE_WIDTH,
    Y_LIM_LOG2,
    REGION_COLORS,
    CNV_COLORS,
    CNV_SHORT_LABELS,
    get_cnv_color,
    format_position,
    sort_chromosomes,
)

logger = logging.getLogger(__name__)


def create_chromosome_plots(
    analysis_result: AnalysisResult,
    ratio_result: Log2RatioResult,
    output_dir: Union[str, Path],
    format: str = "pdf",
    chromosomes: Optional[List[str]] = None,
) -> None:
    """Generate per-chromosome CNV plots.

    Creates two-panel plots for each chromosome:
    - Panel 1: Scatter plot with log2 ratios and segmentation line
    - Panel 2: CNV calls as colored rectangles

    Args:
        analysis_result: Analysis result with segments and CNV calls
        ratio_result: Log2 ratio data with positions
        output_dir: Directory to save plots
        format: Output format ('pdf', 'png', 'svg')
        chromosomes: Specific chromosomes to plot (default: all)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_name = analysis_result.test_sample
    logger.info(f"Creating chromosome plots for {sample_name}")

    # Determine chromosomes to plot
    if chromosomes is None:
        chromosomes = sort_chromosomes(list(analysis_result.chromosome_results.keys()))

    for chrom in chromosomes:
        if chrom not in analysis_result.chromosome_results:
            logger.warning(f"No results for chromosome {chrom}, skipping")
            continue

        chrom_result = analysis_result.chromosome_results[chrom]
        chrom_data = ratio_result.get_chromosome_data(chrom)

        if chrom_data["n_windows"] == 0:
            logger.warning(f"No data for chromosome {chrom}, skipping")
            continue

        output_path = output_dir / f"PlotResults_{chrom}.{format}"
        plot_chromosome(
            chrom=chrom,
            log2_ratios=chrom_data["log2_ratios"],
            positions=chrom_data["positions"],
            in_target_mask=chrom_data["in_target_mask"],
            segments=chrom_result.segments,
            sample_name=sample_name,
            output_path=output_path,
        )
        logger.info(f"  Created PlotResults_{chrom}.{format}")

    logger.info(f"Chromosome plots saved to {output_dir}")


def plot_chromosome(
    chrom: str,
    log2_ratios: np.ndarray,
    positions: np.ndarray,
    in_target_mask: np.ndarray,
    segments: List[CNVSegment],
    sample_name: str = "",
    output_path: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """Create a detailed plot for a single chromosome.

    Creates a two-panel figure:
    - Top: Scatter plot with log2 ratios, colored by IN/OUT target,
           with red segmentation line overlay
    - Bottom: Log2 ratio line with CNV call rectangles

    Args:
        chrom: Chromosome name
        log2_ratios: Log2 ratio values
        positions: Genomic positions (bp)
        in_target_mask: Boolean mask for IN-target windows
        segments: List of CNV segments
        sample_name: Sample name for title
        output_path: Path to save figure (optional)

    Returns:
        Matplotlib Figure object
    """
    fig, axes = plt.subplots(2, 1, figsize=(15, 10), dpi=FIGURE_DPI)

    # Panel 1: Scatter plot with segmentation
    _plot_scatter_panel(
        ax=axes[0],
        positions=positions,
        log2_ratios=log2_ratios,
        in_target_mask=in_target_mask,
        segments=segments,
    )

    # Panel 2: CNV calls
    _plot_cnv_panel(ax=axes[1], positions=positions, log2_ratios=log2_ratios, segments=segments)

    # Title
    fig.suptitle(f"{sample_name} - {chrom}", fontsize=16, fontweight="bold")

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)

    return fig


def _plot_scatter_panel(
    ax: plt.Axes,
    positions: np.ndarray,
    log2_ratios: np.ndarray,
    in_target_mask: np.ndarray,
    segments: List[CNVSegment],
) -> None:
    """Plot scatter panel with log2 ratios and segmentation.

    Args:
        ax: Matplotlib axes
        positions: Genomic positions
        log2_ratios: Log2 ratio values
        in_target_mask: Boolean mask for IN-target
        segments: CNV segments with segment means
    """
    # Separate IN and OUT target points
    out_mask = ~in_target_mask

    # Plot OFF-target first (lighter, in background)
    if np.any(out_mask):
        ax.scatter(
            positions[out_mask],
            log2_ratios[out_mask],
            s=SCATTER_SIZE,
            alpha=SCATTER_ALPHA * 0.7,
            c=REGION_COLORS["OUT"],
            label="Off-target",
            edgecolors="none",
        )

    # Plot IN-target (darker, in foreground)
    if np.any(in_target_mask):
        ax.scatter(
            positions[in_target_mask],
            log2_ratios[in_target_mask],
            s=SCATTER_SIZE,
            alpha=SCATTER_ALPHA,
            c=REGION_COLORS["IN"],
            label="In-target",
            edgecolors="none",
        )

    # Plot segmentation line
    for seg in segments:
        # Draw horizontal line for segment mean
        ax.hlines(
            seg.segment_mean,
            seg.start,
            seg.end,
            colors=SEGMENT_COLOR,
            linewidth=SEGMENT_LINE_WIDTH,
            zorder=10,
        )

    # Reference line at y=0
    ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.7)

    # Styling
    ax.set_xlim(positions.min() - 1e6, positions.max() + 1e6)
    ax.set_ylim(Y_LIM_LOG2)
    ax.set_xlabel("Position (bp)", fontsize=12)
    ax.set_ylabel("log2 ratio", fontsize=12)
    ax.set_title("Segmentation", fontsize=14)

    # Legend
    handles = [
        mpatches.Patch(color=REGION_COLORS["IN"], label="In-target"),
        mpatches.Patch(color=REGION_COLORS["OUT"], label="Off-target"),
        plt.Line2D([0], [0], color=SEGMENT_COLOR, linewidth=2, label="Segment"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=10)

    ax.grid(True, alpha=0.3)

    # Format x-axis with Mb labels
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x/1e6:.0f}"))
    ax.set_xlabel("Position (Mb)", fontsize=12)


def _plot_cnv_panel(
    ax: plt.Axes, positions: np.ndarray, log2_ratios: np.ndarray, segments: List[CNVSegment]
) -> None:
    """Plot CNV calls as colored rectangles.

    Args:
        ax: Matplotlib axes
        positions: Genomic positions
        log2_ratios: Log2 ratio values
        segments: CNV segments with calls
    """
    # Plot log2 ratio as grey line (background)
    ax.plot(positions, log2_ratios, "-", color="grey", linewidth=0.5, alpha=0.5)

    # Reference line at y=0
    ax.axhline(0, color="black", linestyle="-", linewidth=1)

    # Plot CNV call rectangles
    patches = []
    colors = []

    for seg in segments:
        if seg.cn_call == 0:
            # Skip normal segments (don't draw rectangles)
            continue

        # Rectangle from segment start to end
        rect = Rectangle(
            (seg.start, Y_LIM_LOG2[0]),  # bottom-left corner
            seg.end - seg.start,  # width
            Y_LIM_LOG2[1] - Y_LIM_LOG2[0],  # height (full y range)
            alpha=0.4,
        )
        patches.append(rect)
        colors.append(get_cnv_color(seg.cn_call))

    if patches:
        collection = PatchCollection(patches, facecolor=colors, edgecolor="none", alpha=0.4)
        ax.add_collection(collection)

    # Styling
    ax.set_xlim(positions.min() - 1e6, positions.max() + 1e6)
    ax.set_ylim(Y_LIM_LOG2)
    ax.set_xlabel("Position (Mb)", fontsize=12)
    ax.set_ylabel("log2 ratio", fontsize=12)
    ax.set_title("CNV Calls", fontsize=14)

    # Create legend for CNV types (only show types that appear)
    present_calls = set(seg.cn_call for seg in segments if seg.cn_call != 0)
    if present_calls:
        handles = []
        for cn_call in sorted(present_calls):
            color = get_cnv_color(cn_call)
            label = CNV_SHORT_LABELS.get(cn_call, str(cn_call))
            handles.append(mpatches.Patch(color=color, alpha=0.6, label=label))
        ax.legend(handles=handles, loc="upper right", fontsize=10)

    ax.grid(True, alpha=0.3)

    # Format x-axis
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x/1e6:.0f}"))


def plot_chromosome_simple(
    chrom: str,
    log2_ratios: np.ndarray,
    positions: np.ndarray,
    segment_values: np.ndarray,
    output_path: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """Create a simple single-panel chromosome plot.

    Shows log2 ratios with segment overlay.

    Args:
        chrom: Chromosome name
        log2_ratios: Log2 ratio values
        positions: Genomic positions
        segment_values: Segmented values (same length as log2_ratios)
        output_path: Path to save figure (optional)

    Returns:
        Matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=(15, 5), dpi=FIGURE_DPI)

    # Scatter plot
    ax.scatter(
        positions, log2_ratios, s=SCATTER_SIZE, alpha=SCATTER_ALPHA, c="#1f77b4", edgecolors="none"
    )

    # Segment line
    ax.plot(
        positions,
        segment_values,
        "-",
        color=SEGMENT_COLOR,
        linewidth=SEGMENT_LINE_WIDTH,
        label="Segment",
    )

    # Reference line
    ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.7)

    ax.set_xlim(positions.min() - 1e6, positions.max() + 1e6)
    ax.set_ylim(Y_LIM_LOG2)
    ax.set_xlabel("Position (Mb)", fontsize=12)
    ax.set_ylabel("log2 ratio", fontsize=12)
    ax.set_title(f"{chrom}", fontsize=14, fontweight="bold")

    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x/1e6:.0f}"))

    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)

    return fig
