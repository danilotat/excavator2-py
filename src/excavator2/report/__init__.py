"""
EXCAVATOR2 reporting and visualization module.

This module provides functions for generating quality control plots
and CNV visualizations.

QC Plots (from prepare stage):
- Bias plots showing WMRC vs size, mappability, and GC content
- Shows raw vs normalized data to demonstrate correction effectiveness

CNV Plots (from analyze stage):
- Per-chromosome plots with segmentation and CNV calls
- Genome-wide summary plots

Example:
    >>> from excavator2.report import create_qc_plots, create_cnv_plots
    >>> create_qc_plots(sample_data, norm_result, output_dir / "qc")
    >>> create_cnv_plots(analysis_result, output_dir / "plots")
"""

from excavator2.report.qc import (
    create_qc_plots,
    plot_size_bias,
    plot_mappability_bias,
    plot_gc_bias,
)
from excavator2.report.chromosome import (
    create_chromosome_plots,
    plot_chromosome,
)
from excavator2.report.genome import (
    create_genome_plot,
    create_genome_heatmap,
    create_cnv_summary_table,
)

__all__ = [
    # QC plots
    "create_qc_plots",
    "plot_size_bias",
    "plot_mappability_bias",
    "plot_gc_bias",
    # Chromosome plots
    "create_chromosome_plots",
    "plot_chromosome",
    # Genome-wide plots
    "create_genome_plot",
    "create_genome_heatmap",
    "create_cnv_summary_table",
]
