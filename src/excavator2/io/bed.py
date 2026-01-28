"""
BED file and coordinate file parsing utilities.

This module provides functions for reading BED files, chromosome coordinate
files, and gap/centromere annotation files used in the EXCAVATOR2 pipeline.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Optional, Union, Tuple
import logging

logger = logging.getLogger(__name__)


@dataclass
class TargetRegion:
    """A target region from a BED file.

    Attributes:
        chrom: Chromosome name
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        name: Optional region name (e.g., gene name)
        score: Optional score field
        strand: Optional strand (+/-)
    """
    chrom: str
    start: int
    end: int
    name: str = ""
    score: float = 0.0
    strand: str = "."

    @property
    def length(self) -> int:
        """Region length in base pairs."""
        return self.end - self.start

    @property
    def midpoint(self) -> int:
        """Midpoint position (0-based)."""
        return (self.start + self.end) // 2


@dataclass
class ChromosomeInfo:
    """Chromosome size information.

    Attributes:
        name: Chromosome name (e.g., "chr1" or "1")
        start: Start position (usually 0 or 1)
        end: End position (chromosome length)
    """
    name: str
    start: int
    end: int

    @property
    def length(self) -> int:
        """Chromosome length."""
        return self.end - self.start


@dataclass
class GapRegion:
    """A genomic gap region (centromere, telomere, etc.).

    Attributes:
        chrom: Chromosome name
        start: Start position (0-based)
        end: End position (0-based, exclusive)
        name: Gap type (e.g., "centromere", "telomere")
    """
    chrom: str
    start: int
    end: int
    name: str = ""

    @property
    def length(self) -> int:
        """Gap length in base pairs."""
        return self.end - self.start


def load_target_bed(
    bed_path: Union[str, Path],
    skip_header: bool = False,
    min_length: int = 0,
    chromosomes: Optional[List[str]] = None
) -> List[TargetRegion]:
    """Load target regions from a BED file.

    Args:
        bed_path: Path to BED file
        skip_header: Skip first line if it's a header
        min_length: Minimum region length to include
        chromosomes: If provided, only include these chromosomes

    Returns:
        List of TargetRegion objects sorted by chromosome and position

    Note:
        BED format is 0-based, half-open [start, end)
        Supports BED3 through BED6 formats
    """
    regions = []
    bed_path = Path(bed_path)

    # Build chromosome set for filtering (handle both chr-prefixed and non-prefixed)
    chrom_set = None
    if chromosomes:
        chrom_set = set(chromosomes)
        # Also add alternate naming conventions
        for chrom in chromosomes:
            if chrom.startswith('chr'):
                chrom_set.add(chrom[3:])
            else:
                chrom_set.add(f'chr{chrom}')

    with open(bed_path) as f:
        for i, line in enumerate(f):
            if skip_header and i == 0:
                continue

            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue

            fields = line.split('\t')
            if len(fields) < 3:
                continue

            chrom = fields[0]

            # Filter by chromosome if specified
            if chrom_set is not None:
                if chrom not in chrom_set:
                    continue

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                logger.warning(f"Skipping invalid line {i+1}: cannot parse coordinates")
                continue

            # Filter by minimum length
            if end - start < min_length:
                continue

            name = fields[3] if len(fields) > 3 else ""
            score = float(fields[4]) if len(fields) > 4 else 0.0
            strand = fields[5] if len(fields) > 5 else "."

            regions.append(TargetRegion(chrom, start, end, name, score, strand))

    # Sort by chromosome and position
    regions.sort(key=lambda r: (r.chrom, r.start, r.end))

    logger.info(f"Loaded {len(regions)} target regions from {bed_path}")
    return regions


def load_chromosome_coordinates(
    coord_path: Union[str, Path],
    skip_header: bool = False
) -> Dict[str, ChromosomeInfo]:
    """Load chromosome coordinate information.

    Args:
        coord_path: Path to chromosome coordinates file
        skip_header: Skip first line (default: False)

    Returns:
        Dictionary mapping chromosome name to ChromosomeInfo

    Expected format (BED):
        chr1  1      248956422
        chr2  1      242193529
        ...
    """
    chromosomes = {}
    coord_path = Path(coord_path)

    with open(coord_path) as f:
        for i, line in enumerate(f):
            if skip_header and i == 0:
                continue

            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 3:
                continue

            name = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                logger.warning(f"Skipping invalid line {i+1}: cannot parse coordinates")
                continue

            chromosomes[name] = ChromosomeInfo(name, start, end)

    logger.info(f"Loaded {len(chromosomes)} chromosome coordinates from {coord_path}")
    return chromosomes


def load_gap_file(
    gap_path: Union[str, Path],
    skip_header: bool = True,
    exclude_alt: bool = True
) -> List[GapRegion]:
    """Load genomic gap regions (centromeres, telomeres, etc.).

    Args:
        gap_path: Path to gap file (UCSC format)
        skip_header: Skip first line (default: True)
        exclude_alt: Exclude alt/random chromosomes containing "_" (default: True)

    Returns:
        List of GapRegion objects

    Expected UCSC gap format (tab-delimited):
        bin chrom chromStart chromEnd ix n size type bridge
        585 chr1  0          10000    1  N 10000 telomere no
        ...

    We use columns: chrom (1), chromStart (2), chromEnd (3), type (7)
    """
    gaps = []
    gap_path = Path(gap_path)

    with open(gap_path) as f:
        for i, line in enumerate(f):
            if skip_header and i == 0:
                continue

            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 4:
                continue

            chrom = fields[1] if len(fields) > 1 else fields[0]

            # Exclude alt/random chromosomes (contain underscore)
            if exclude_alt and '_' in chrom:
                continue

            try:
                start = int(fields[2])
                end = int(fields[3])
            except ValueError:
                logger.warning(f"Skipping invalid gap line {i+1}")
                continue

            # Gap type is in column 7 if present
            name = fields[7] if len(fields) > 7 else ""

            gaps.append(GapRegion(chrom, start, end, name))

    # Sort by chromosome and position
    gaps.sort(key=lambda g: (g.chrom, g.start))

    logger.info(f"Loaded {len(gaps)} gap regions from {gap_path}")
    return gaps


def regions_to_bed(
    regions: List[TargetRegion],
    output_path: Union[str, Path],
    include_header: bool = False
) -> None:
    """Write regions to a BED file.

    Args:
        regions: List of TargetRegion objects
        output_path: Output file path
        include_header: Include header line
    """
    output_path = Path(output_path)

    with open(output_path, 'w') as f:
        if include_header:
            f.write("#chrom\tstart\tend\tname\tscore\tstrand\n")

        for region in regions:
            f.write(f"{region.chrom}\t{region.start}\t{region.end}")
            if region.name or region.score or region.strand != ".":
                f.write(f"\t{region.name}\t{region.score}\t{region.strand}")
            f.write("\n")

    logger.info(f"Wrote {len(regions)} regions to {output_path}")


def merge_overlapping_regions(
    regions: List[TargetRegion],
    merge_distance: int = 0
) -> List[TargetRegion]:
    """Merge overlapping or adjacent target regions.

    Args:
        regions: List of TargetRegion objects (will be sorted)
        merge_distance: Merge regions within this distance (default: 0 = touching)

    Returns:
        List of merged TargetRegion objects
    """
    if not regions:
        return []

    # Sort by chromosome and position
    sorted_regions = sorted(regions, key=lambda r: (r.chrom, r.start, r.end))

    merged = []
    current = sorted_regions[0]

    for region in sorted_regions[1:]:
        # Same chromosome and overlapping/adjacent
        if region.chrom == current.chrom and region.start <= current.end + merge_distance:
            # Extend current region
            current = TargetRegion(
                chrom=current.chrom,
                start=current.start,
                end=max(current.end, region.end),
                name=current.name or region.name
            )
        else:
            merged.append(current)
            current = region

    merged.append(current)

    logger.info(f"Merged {len(regions)} regions into {len(merged)} non-overlapping regions")
    return merged


def get_chromosome_regions(
    regions: List[TargetRegion],
    chromosome: str
) -> List[TargetRegion]:
    """Get regions for a specific chromosome.

    Args:
        regions: List of TargetRegion objects
        chromosome: Chromosome name to filter by

    Returns:
        List of TargetRegion objects for the specified chromosome
    """
    # Handle both chr-prefixed and non-prefixed naming
    chrom_variants = {chromosome}
    if chromosome.startswith('chr'):
        chrom_variants.add(chromosome[3:])
    else:
        chrom_variants.add(f'chr{chromosome}')

    return [r for r in regions if r.chrom in chrom_variants]


def normalize_chromosome_name(
    chrom: str,
    use_chr_prefix: bool = True
) -> str:
    """Normalize chromosome name to consistent format.

    Args:
        chrom: Chromosome name (with or without 'chr' prefix)
        use_chr_prefix: If True, add 'chr' prefix; if False, remove it

    Returns:
        Normalized chromosome name
    """
    if use_chr_prefix:
        return chrom if chrom.startswith('chr') else f'chr{chrom}'
    else:
        return chrom[3:] if chrom.startswith('chr') else chrom
