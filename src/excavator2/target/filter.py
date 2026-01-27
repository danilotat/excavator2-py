"""
Target region filtering and window creation.

This module provides functions for:
- Creating analysis windows from target BED regions
- Filtering windows that overlap genomic gaps (centromeres, telomeres)
- Creating both IN-target (exon) and OUT-target (inter-exon) windows
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Optional, Union, Tuple
import logging

import numpy as np

from excavator2.io.bed import (
    TargetRegion,
    ChromosomeInfo,
    GapRegion,
    load_target_bed,
    load_gap_file,
    load_chromosome_coordinates,
)

logger = logging.getLogger(__name__)


@dataclass
class AnalysisWindow:
    """A window for CNV analysis.

    Attributes:
        chrom: Chromosome name
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        window_id: Unique window identifier
        region_class: 'IN' for target regions, 'OUT' for inter-target
        position: Center position of the window
    """
    chrom: str
    start: int
    end: int
    window_id: str
    region_class: str  # 'IN' or 'OUT'
    position: int = 0

    def __post_init__(self):
        if self.position == 0:
            self.position = (self.start + self.end) // 2

    @property
    def length(self) -> int:
        """Window length in base pairs."""
        return self.end - self.start


@dataclass
class FilteredTargetResult:
    """Result of target filtering.

    Attributes:
        windows: List of analysis windows
        chromosomes: List of unique chromosomes
        n_in_target: Number of IN-target windows
        n_out_target: Number of OUT-target windows
        n_filtered: Number of windows filtered out (gap overlap)
        target_name: Name of the target
    """
    windows: List[AnalysisWindow]
    chromosomes: List[str]
    n_in_target: int
    n_out_target: int
    n_filtered: int
    target_name: str

    @property
    def n_windows(self) -> int:
        """Total number of windows."""
        return len(self.windows)

    def get_chromosome_windows(self, chromosome: str) -> List[AnalysisWindow]:
        """Get windows for a specific chromosome."""
        # Handle both chr-prefixed and non-prefixed naming
        chrom_variants = {chromosome}
        if chromosome.startswith('chr'):
            chrom_variants.add(chromosome[3:])
        else:
            chrom_variants.add(f'chr{chromosome}')

        return [w for w in self.windows if w.chrom in chrom_variants]


def create_analysis_windows(
    target_regions: List[TargetRegion],
    chromosome_info: Dict[str, ChromosomeInfo],
    window_size: int,
    flank: int = 200
) -> List[AnalysisWindow]:
    """Create analysis windows from target regions.

    This function creates both IN-target windows (for exons) and OUT-target
    windows (for inter-exon gaps). OUT-target windows are only created when
    the gap between exons is large enough.

    Args:
        target_regions: List of target regions (exons)
        chromosome_info: Chromosome size information
        window_size: Size of each analysis window
        flank: Buffer size around gaps (default: 200bp)

    Returns:
        List of AnalysisWindow objects
    """
    windows = []
    window_counter = 1

    # Group regions by chromosome
    regions_by_chrom: Dict[str, List[TargetRegion]] = {}
    for region in target_regions:
        if region.chrom not in regions_by_chrom:
            regions_by_chrom[region.chrom] = []
        regions_by_chrom[region.chrom].append(region)

    for chrom in sorted(regions_by_chrom.keys()):
        # Skip chromosomes not in chromosome_info
        if chrom not in chromosome_info:
            # Try alternate naming
            alt_chrom = f'chr{chrom}' if not chrom.startswith('chr') else chrom[3:]
            if alt_chrom not in chromosome_info:
                logger.warning(f"Chromosome {chrom} not in chromosome info, skipping")
                continue

        chrom_regions = sorted(regions_by_chrom[chrom], key=lambda r: r.start)

        # Process each region and gaps between them
        prev_end = 0

        for region in chrom_regions:
            # Create OUT-target windows for gap before this region
            gap_size = region.start - prev_end
            if gap_size >= window_size + 2 * flank:
                # Calculate number of windows that fit
                usable_gap = gap_size - 2 * flank
                n_out_windows = usable_gap // window_size

                for j in range(n_out_windows):
                    out_start = prev_end + flank + j * window_size
                    out_end = out_start + window_size

                    windows.append(AnalysisWindow(
                        chrom=chrom,
                        start=out_start,
                        end=out_end,
                        window_id=f"a{window_counter}",
                        region_class="OUT"
                    ))
                    window_counter += 1

            # Create IN-target window for the exon
            windows.append(AnalysisWindow(
                chrom=chrom,
                start=region.start,
                end=region.end,
                window_id=f"a{window_counter}",
                region_class="IN"
            ))
            window_counter += 1

            prev_end = region.end

    logger.info(f"Created {len(windows)} analysis windows")
    return windows


def filter_gap_overlapping_windows(
    windows: List[AnalysisWindow],
    gaps: List[GapRegion]
) -> Tuple[List[AnalysisWindow], int]:
    """Filter windows that overlap genomic gaps.

    Args:
        windows: List of analysis windows
        gaps: List of gap regions (centromeres, telomeres, etc.)

    Returns:
        Tuple of (filtered_windows, n_filtered)
    """
    if not gaps:
        return windows, 0

    # Build gap lookup by chromosome
    gaps_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    for gap in gaps:
        if gap.chrom not in gaps_by_chrom:
            gaps_by_chrom[gap.chrom] = []
        gaps_by_chrom[gap.chrom].append((gap.start, gap.end))

    # Sort gaps by start position
    for chrom in gaps_by_chrom:
        gaps_by_chrom[chrom].sort()

    filtered_windows = []
    n_filtered = 0

    for window in windows:
        # Get gaps for this chromosome
        chrom_gaps = gaps_by_chrom.get(window.chrom, [])
        if not chrom_gaps:
            # Try alternate naming
            alt_chrom = f'chr{window.chrom}' if not window.chrom.startswith('chr') else window.chrom[3:]
            chrom_gaps = gaps_by_chrom.get(alt_chrom, [])

        # Check overlap with any gap
        overlaps_gap = False
        for gap_start, gap_end in chrom_gaps:
            # Check if window overlaps gap
            if window.start < gap_end and window.end > gap_start:
                overlaps_gap = True
                break

        if overlaps_gap:
            n_filtered += 1
        else:
            filtered_windows.append(window)

    logger.info(f"Filtered {n_filtered} windows overlapping gaps, {len(filtered_windows)} remaining")
    return filtered_windows, n_filtered


def create_filtered_target(
    bed_path: Union[str, Path],
    chromosome_path: Union[str, Path],
    gap_path: Union[str, Path],
    window_size: int,
    target_name: str,
    flank: int = 200,
    min_region_length: int = 0,
    chromosomes: Optional[List[str]] = None
) -> FilteredTargetResult:
    """Create filtered analysis windows from target BED file.

    This is the main entry point for target preparation. It:
    1. Loads target regions from BED file
    2. Loads chromosome coordinates
    3. Creates analysis windows (IN-target and OUT-target)
    4. Filters windows overlapping genomic gaps

    Args:
        bed_path: Path to target BED file
        chromosome_path: Path to chromosome coordinates file
        gap_path: Path to gap annotations file
        window_size: Size of analysis windows
        target_name: Name for this target set
        flank: Buffer around gaps (default: 200bp)
        min_region_length: Minimum region length to include
        chromosomes: If provided, only process these chromosomes

    Returns:
        FilteredTargetResult with analysis windows
    """
    logger.info(f"Creating filtered target from {bed_path}")
    logger.info(f"Window size: {window_size}, Flank: {flank}")

    # Load target regions
    target_regions = load_target_bed(
        bed_path,
        min_length=min_region_length,
        chromosomes=chromosomes
    )
    logger.info(f"Loaded {len(target_regions)} target regions")

    # Load chromosome coordinates
    chromosome_info = load_chromosome_coordinates(chromosome_path)
    logger.info(f"Loaded {len(chromosome_info)} chromosome coordinates")

    # Load gap regions
    gaps = load_gap_file(gap_path)
    logger.info(f"Loaded {len(gaps)} gap regions")

    # Create analysis windows
    windows = create_analysis_windows(
        target_regions,
        chromosome_info,
        window_size,
        flank
    )

    # Count IN/OUT before filtering
    n_in_before = sum(1 for w in windows if w.region_class == "IN")
    n_out_before = sum(1 for w in windows if w.region_class == "OUT")

    # Filter gap-overlapping windows
    filtered_windows, n_filtered = filter_gap_overlapping_windows(windows, gaps)

    # Count IN/OUT after filtering
    n_in_target = sum(1 for w in filtered_windows if w.region_class == "IN")
    n_out_target = sum(1 for w in filtered_windows if w.region_class == "OUT")

    # Get unique chromosomes
    unique_chroms = sorted(set(w.chrom for w in filtered_windows))

    logger.info(f"Final windows: {len(filtered_windows)} ({n_in_target} IN, {n_out_target} OUT)")
    logger.info(f"Filtered {n_filtered} windows overlapping gaps")

    return FilteredTargetResult(
        windows=filtered_windows,
        chromosomes=unique_chroms,
        n_in_target=n_in_target,
        n_out_target=n_out_target,
        n_filtered=n_filtered,
        target_name=target_name
    )


def renumber_windows(windows: List[AnalysisWindow]) -> List[AnalysisWindow]:
    """Renumber window IDs sequentially.

    Args:
        windows: List of analysis windows

    Returns:
        List of windows with sequential IDs (a1, a2, a3, ...)
    """
    renumbered = []
    for i, window in enumerate(windows, start=1):
        renumbered.append(AnalysisWindow(
            chrom=window.chrom,
            start=window.start,
            end=window.end,
            window_id=f"a{i}",
            region_class=window.region_class,
            position=window.position
        ))
    return renumbered
