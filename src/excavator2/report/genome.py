"""
Genome-wide CNV visualization.

This module generates genome-wide plots showing CNV calls across
all chromosomes in a single view. Uses a linear layout with
chromosome ideograms.
"""

from pathlib import Path
from typing import List, Optional, Union, Dict
import logging

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.collections import PatchCollection

from excavator2.analyze.pipeline import AnalysisResult, CNVSegment
from excavator2.report.utils import (
    FIGURE_DPI,
    CNV_COLORS,
    CNV_SHORT_LABELS,
    CHROMOSOME_ORDER,
    CHROMOSOME_LENGTHS_HG38,
    get_cnv_color,
    sort_chromosomes,
    normalize_chromosome_name,
)

logger = logging.getLogger(__name__)


def create_genome_plot(
    analysis_result: AnalysisResult,
    output_path: Optional[Union[str, Path]] = None,
    assembly: str = "hg38",
    show_normal: bool = False,
    title: Optional[str] = None,
) -> plt.Figure:
    """Create a genome-wide CNV plot.

    Shows all chromosomes in a linear layout with CNV calls
    displayed as colored bars.

    Args:
        analysis_result: Analysis result with CNV segments
        output_path: Path to save figure (optional)
        assembly: Genome assembly ('hg38', 'hg19')
        show_normal: Show normal segments (default: False)
        title: Custom title (default: sample name)

    Returns:
        Matplotlib Figure object
    """
    sample_name = analysis_result.test_sample

    if title is None:
        title = f"EXCAVATOR2 Genome-wide CNV Calls - {sample_name}"

    logger.info(f"Creating genome-wide plot for {sample_name}")

    # Get chromosome lengths
    chrom_lengths = CHROMOSOME_LENGTHS_HG38.copy()

    # Determine which chromosomes have data
    available_chroms = []
    for chrom in CHROMOSOME_ORDER:
        norm_chrom = normalize_chromosome_name(chrom)
        # Check both prefixed and unprefixed versions
        if chrom in analysis_result.chromosome_results:
            available_chroms.append(chrom)
        elif chrom.replace("chr", "") in analysis_result.chromosome_results:
            available_chroms.append(chrom.replace("chr", ""))
        elif norm_chrom in analysis_result.chromosome_results:
            available_chroms.append(norm_chrom)

    if not available_chroms:
        logger.warning("No chromosome data available for genome-wide plot")
        fig, ax = plt.subplots(figsize=(15, 4))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=14)
        if output_path:
            fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
        return fig

    # Create figure
    fig, ax = plt.subplots(figsize=(18, 8), dpi=FIGURE_DPI)

    # Calculate genome positions (cumulative)
    genome_offset = {}
    current_offset = 0
    for chrom in CHROMOSOME_ORDER:
        norm_chrom = normalize_chromosome_name(chrom)
        if norm_chrom in chrom_lengths:
            genome_offset[chrom] = current_offset
            genome_offset[norm_chrom] = current_offset
            genome_offset[chrom.replace("chr", "")] = current_offset
            current_offset += chrom_lengths[norm_chrom]

    total_genome_length = current_offset

    # Draw chromosome ideograms
    ideogram_height = 0.3
    ideogram_y = 0.5

    for i, chrom in enumerate(CHROMOSOME_ORDER):
        norm_chrom = normalize_chromosome_name(chrom)
        if norm_chrom not in chrom_lengths:
            continue

        chrom_start = genome_offset.get(chrom, 0)
        chrom_len = chrom_lengths[norm_chrom]

        # Alternating colors for chromosomes
        color = "#e0e0e0" if i % 2 == 0 else "#c0c0c0"

        rect = FancyBboxPatch(
            (chrom_start, ideogram_y - ideogram_height / 2),
            chrom_len,
            ideogram_height,
            boxstyle="round,pad=0,rounding_size=0",
            facecolor=color,
            edgecolor="black",
            linewidth=0.5,
        )
        ax.add_patch(rect)

        # Chromosome label
        ax.text(
            chrom_start + chrom_len / 2,
            ideogram_y - ideogram_height / 2 - 0.15,
            chrom.replace("chr", ""),
            ha="center",
            va="top",
            fontsize=8,
            rotation=45,
        )

    # Draw CNV calls
    cnv_bar_height = 0.15
    del_y = ideogram_y + ideogram_height / 2 + 0.1  # Above ideogram
    amp_y = ideogram_y - ideogram_height / 2 - 0.1 - cnv_bar_height  # Below ideogram

    del_patches = []
    del_colors = []
    amp_patches = []
    amp_colors = []

    for chrom, chrom_result in analysis_result.chromosome_results.items():
        # Find offset for this chromosome
        offset = None
        for key in [chrom, normalize_chromosome_name(chrom), f"chr{chrom}"]:
            if key in genome_offset:
                offset = genome_offset[key]
                break

        if offset is None:
            continue

        for seg in chrom_result.segments:
            if seg.cn_call == 0 and not show_normal:
                continue

            genome_start = offset + seg.start
            width = seg.end - seg.start

            if seg.cn_call < 0:  # Deletion
                y_pos = del_y + (0.2 if seg.cn_call == -2 else 0)
                rect = Rectangle((genome_start, y_pos), width, cnv_bar_height)
                del_patches.append(rect)
                del_colors.append(get_cnv_color(seg.cn_call))
            elif seg.cn_call > 0:  # Gain
                y_pos = amp_y - (0.2 if seg.cn_call == 2 else 0)
                rect = Rectangle((genome_start, y_pos), width, cnv_bar_height)
                amp_patches.append(rect)
                amp_colors.append(get_cnv_color(seg.cn_call))

    # Add collections
    if del_patches:
        del_collection = PatchCollection(
            del_patches, facecolor=del_colors, edgecolor="none", alpha=0.8
        )
        ax.add_collection(del_collection)

    if amp_patches:
        amp_collection = PatchCollection(
            amp_patches, facecolor=amp_colors, edgecolor="none", alpha=0.8
        )
        ax.add_collection(amp_collection)

    # Styling
    ax.set_xlim(-total_genome_length * 0.02, total_genome_length * 1.02)
    ax.set_ylim(-0.8, 1.2)

    ax.set_xlabel("Genome Position", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=20)

    # Hide y-axis
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Format x-axis (Gb)
    def format_gb(x, p):
        return f"{x / 1e9:.1f} Gb"

    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_gb))

    # Add labels for gain/loss regions
    ax.text(
        -total_genome_length * 0.015,
        del_y + cnv_bar_height / 2,
        "DEL",
        ha="right",
        va="center",
        fontsize=10,
        fontweight="bold",
        color=CNV_COLORS[-1],
    )
    ax.text(
        -total_genome_length * 0.015,
        amp_y + cnv_bar_height / 2,
        "AMP",
        ha="right",
        va="center",
        fontsize=10,
        fontweight="bold",
        color=CNV_COLORS[1],
    )

    # Legend
    handles = []
    for cn_call in [-2, -1, 1, 2]:
        color = get_cnv_color(cn_call)
        label = CNV_SHORT_LABELS.get(cn_call, str(cn_call))
        handles.append(mpatches.Patch(color=color, alpha=0.8, label=label))

    ax.legend(handles=handles, loc="upper right", ncol=4, fontsize=10, framealpha=0.9)

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"Genome-wide plot saved to {output_path}")

    return fig


def create_genome_heatmap(
    analysis_results: List[AnalysisResult],
    output_path: Optional[Union[str, Path]] = None,
    title: str = "CNV Heatmap",
) -> plt.Figure:
    """Create a heatmap showing CNVs across multiple samples.

    Each row is a sample, each column represents a genomic bin.
    Colors indicate copy number state.

    Args:
        analysis_results: List of analysis results for multiple samples
        output_path: Path to save figure (optional)
        title: Plot title

    Returns:
        Matplotlib Figure object
    """
    if not analysis_results:
        fig, ax = plt.subplots(figsize=(15, 4))
        ax.text(0.5, 0.5, "No data available", ha="center", va="center")
        return fig

    logger.info(f"Creating CNV heatmap for {len(analysis_results)} samples")

    # Determine bin size for genome discretization
    bin_size = 1_000_000  # 1 Mb bins
    total_bins = sum(length // bin_size + 1 for length in CHROMOSOME_LENGTHS_HG38.values())

    # Create chromosome offset map
    chrom_offsets = {}
    current_offset = 0
    for chrom in CHROMOSOME_ORDER:
        norm_chrom = normalize_chromosome_name(chrom)
        if norm_chrom in CHROMOSOME_LENGTHS_HG38:
            chrom_offsets[chrom] = current_offset
            chrom_offsets[norm_chrom] = current_offset
            n_bins = CHROMOSOME_LENGTHS_HG38[norm_chrom] // bin_size + 1
            current_offset += n_bins

    # Create data matrix
    n_samples = len(analysis_results)
    data = np.zeros((n_samples, current_offset))
    sample_names = []

    for i, result in enumerate(analysis_results):
        sample_names.append(result.test_sample)

        for chrom, chrom_result in result.chromosome_results.items():
            # Get offset
            offset = None
            for key in [chrom, normalize_chromosome_name(chrom)]:
                if key in chrom_offsets:
                    offset = chrom_offsets[key]
                    break
            if offset is None:
                continue

            for seg in chrom_result.segments:
                start_bin = offset + seg.start // bin_size
                end_bin = offset + seg.end // bin_size
                data[i, start_bin : end_bin + 1] = seg.cn_call

    # Create figure
    fig, ax = plt.subplots(figsize=(18, max(4, n_samples * 0.5)), dpi=FIGURE_DPI)

    # Custom colormap for CNV states
    from matplotlib.colors import ListedColormap, BoundaryNorm

    colors_list = [
        CNV_COLORS[-2],  # 2-copy DEL
        CNV_COLORS[-1],  # 1-copy DEL
        (0.95, 0.95, 0.95),  # Normal
        CNV_COLORS[1],  # 1-copy AMP
        CNV_COLORS[2],  # N-copy AMP
    ]
    cmap = ListedColormap(colors_list)
    bounds = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
    norm = BoundaryNorm(bounds, cmap.N)

    # Plot heatmap
    im = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")

    # Y-axis labels (sample names)
    ax.set_yticks(range(n_samples))
    ax.set_yticklabels(sample_names, fontsize=10)

    # Add chromosome boundaries
    prev_offset = 0
    for chrom in CHROMOSOME_ORDER:
        if chrom in chrom_offsets:
            norm_chrom = normalize_chromosome_name(chrom)
            n_bins = CHROMOSOME_LENGTHS_HG38.get(norm_chrom, 0) // bin_size + 1
            offset = chrom_offsets[chrom] + n_bins
            ax.axvline(offset, color="grey", linewidth=0.5, alpha=0.5)

            # Chromosome label
            mid_x = prev_offset + n_bins / 2
            ax.text(
                mid_x,
                -0.7,
                chrom.replace("chr", ""),
                ha="center",
                va="top",
                fontsize=7,
                rotation=45,
            )
            prev_offset = offset

    ax.set_xlabel("Genome Position", fontsize=12)
    ax.set_ylabel("Sample", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_ticks([-2, -1, 0, 1, 2])
    cbar.set_ticklabels(["2-DEL", "DEL", "Normal", "AMP", "2-AMP"])

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"CNV heatmap saved to {output_path}")

    return fig


def create_cnv_summary_table(
    analysis_result: AnalysisResult, output_path: Optional[Union[str, Path]] = None
) -> plt.Figure:
    """Create a visual summary table of CNV calls.

    Shows a table with CNV statistics per chromosome.

    Args:
        analysis_result: Analysis result
        output_path: Path to save figure (optional)

    Returns:
        Matplotlib Figure object
    """
    sample_name = analysis_result.test_sample

    # Collect statistics per chromosome
    stats = []
    for chrom in sort_chromosomes(list(analysis_result.chromosome_results.keys())):
        chrom_result = analysis_result.chromosome_results[chrom]

        n_segments = chrom_result.n_segments
        n_cnvs = chrom_result.n_cnvs
        n_del = sum(1 for s in chrom_result.segments if s.cn_call < 0)
        n_amp = sum(1 for s in chrom_result.segments if s.cn_call > 0)

        # Calculate total CNV length
        del_length = sum(s.length for s in chrom_result.segments if s.cn_call < 0)
        amp_length = sum(s.length for s in chrom_result.segments if s.cn_call > 0)

        stats.append(
            {
                "Chromosome": chrom,
                "Segments": n_segments,
                "CNVs": n_cnvs,
                "Deletions": n_del,
                "Del Length (Mb)": del_length / 1e6,
                "Gains": n_amp,
                "Gain Length (Mb)": amp_length / 1e6,
            }
        )

    # Create figure with table
    fig, ax = plt.subplots(figsize=(12, max(4, len(stats) * 0.3 + 2)), dpi=FIGURE_DPI)
    ax.axis("off")

    # Create table
    columns = ["Chromosome", "Segments", "CNVs", "Deletions", "Del (Mb)", "Gains", "Gain (Mb)"]
    cell_text = []
    for s in stats:
        cell_text.append(
            [
                s["Chromosome"],
                str(s["Segments"]),
                str(s["CNVs"]),
                str(s["Deletions"]),
                f"{s['Del Length (Mb)']:.2f}",
                str(s["Gains"]),
                f"{s['Gain Length (Mb)']:.2f}",
            ]
        )

    # Add totals row
    total_segs = sum(s["Segments"] for s in stats)
    total_cnvs = sum(s["CNVs"] for s in stats)
    total_del = sum(s["Deletions"] for s in stats)
    total_del_len = sum(s["Del Length (Mb)"] for s in stats)
    total_amp = sum(s["Gains"] for s in stats)
    total_amp_len = sum(s["Gain Length (Mb)"] for s in stats)

    cell_text.append(
        [
            "TOTAL",
            str(total_segs),
            str(total_cnvs),
            str(total_del),
            f"{total_del_len:.2f}",
            str(total_amp),
            f"{total_amp_len:.2f}",
        ]
    )

    table = ax.table(cellText=cell_text, colLabels=columns, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)

    # Style header
    for j, col in enumerate(columns):
        table[(0, j)].set_facecolor("#4a7ebb")
        table[(0, j)].set_text_props(color="white", fontweight="bold")

    # Style totals row
    for j in range(len(columns)):
        table[(len(cell_text), j)].set_facecolor("#e0e0e0")
        table[(len(cell_text), j)].set_text_props(fontweight="bold")

    ax.set_title(f"CNV Summary - {sample_name}", fontsize=14, fontweight="bold", pad=20)

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)

    return fig
