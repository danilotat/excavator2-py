"""
EXCAVATOR2 target initialization module.

This module provides functionality for:
- Creating analysis windows from target BED regions
- Filtering windows overlapping genomic gaps (centromeres, telomeres)
- Calculating GC content from reference genome
- Extracting mappability from BigWig files
- Saving/loading target data in HDF5 format

Example:
    >>> from excavator2.target import initialize_target, save_target_data
    >>> target = initialize_target(
    ...     bed_path="targets.bed",
    ...     fasta_path="ref.fasta",
    ...     bigwig_path="mappability.bw",
    ...     chromosome_path="chromosomes.txt",
    ...     gap_path="gaps.txt",
    ...     target_name="MyTarget",
    ...     assembly="hg38",
    ...     window_size=30000
    ... )
    >>> save_target_data(target, "target.h5")
"""

# Window filtering
from excavator2.target.filter import (
    AnalysisWindow,
    FilteredTargetResult,
    create_analysis_windows,
    filter_gap_overlapping_windows,
    create_filtered_target,
    renumber_windows,
)

# GC content
from excavator2.target.gc_content import (
    calculate_gc_content,
    calculate_gc_content_by_chromosome,
    get_window_gc_stats,
)

# Mappability
from excavator2.target.mappability import (
    extract_mappability,
    extract_mappability_by_chromosome,
    get_mappability_stats,
    filter_low_mappability_windows,
)

# Main initialization
from excavator2.target.init import (
    TargetData,
    initialize_target,
    save_target_data,
    load_target_data,
)

__all__ = [
    # Window filtering
    "AnalysisWindow",
    "FilteredTargetResult",
    "create_analysis_windows",
    "filter_gap_overlapping_windows",
    "create_filtered_target",
    "renumber_windows",
    # GC content
    "calculate_gc_content",
    "calculate_gc_content_by_chromosome",
    "get_window_gc_stats",
    # Mappability
    "extract_mappability",
    "extract_mappability_by_chromosome",
    "get_mappability_stats",
    "filter_low_mappability_windows",
    # Main initialization
    "TargetData",
    "initialize_target",
    "save_target_data",
    "load_target_data",
]
