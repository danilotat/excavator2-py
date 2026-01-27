"""
BAM/CRAM file I/O utilities.

This module provides functions for reading BAM/CRAM files and counting
reads in genomic regions using pysam.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union, Iterator, Tuple
import numpy as np
import pysam


@dataclass
class GenomicRegion:
    """A genomic region (window/exon).

    Attributes:
        chrom: Chromosome name
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        name: Optional region name (e.g., gene name)
    """
    chrom: str
    start: int
    end: int
    name: str = ""

    @property
    def length(self) -> int:
        """Region length in base pairs."""
        return self.end - self.start

    @property
    def midpoint(self) -> int:
        """Midpoint position."""
        return (self.start + self.end) // 2


@dataclass
class ReadCountResult:
    """Result of read counting for a set of regions.

    Attributes:
        counts: Read counts per region
        regions: List of regions that were counted
        total_reads: Total reads processed
        filtered_reads: Reads filtered out (low MAPQ, duplicates, etc.)
    """
    counts: np.ndarray
    regions: List[GenomicRegion]
    total_reads: int
    filtered_reads: int

    @property
    def n_regions(self) -> int:
        """Number of regions."""
        return len(self.regions)

    @property
    def mean_count(self) -> float:
        """Mean read count per region."""
        if len(self.counts) == 0:
            return 0.0
        return float(np.mean(self.counts))


class BAMReader:
    """BAM/CRAM file reader for read counting.

    This class provides methods to count reads in genomic regions from
    BAM or CRAM files using pysam. It filters reads by mapping quality
    and excludes unmapped reads and duplicates.

    Example:
        >>> reader = BAMReader("sample.bam", min_mapq=20)
        >>> regions = [GenomicRegion("chr1", 1000, 2000)]
        >>> result = reader.count_reads(regions)
        >>> print(f"Counted {result.counts[0]} reads in region")

    Args:
        bam_path: Path to BAM or CRAM file (must be indexed)
        min_mapq: Minimum mapping quality (default: 20)
        reference: Path to reference FASTA (required for CRAM)
    """

    # SAM flags
    FLAG_UNMAPPED = 0x4
    FLAG_DUPLICATE = 0x400
    EXCLUDE_FLAGS = FLAG_UNMAPPED | FLAG_DUPLICATE  # 1028

    def __init__(
        self,
        bam_path: Union[str, Path],
        min_mapq: int = 20,
        reference: Optional[Union[str, Path]] = None
    ):
        self.bam_path = Path(bam_path)
        self.min_mapq = min_mapq
        self.reference = Path(reference) if reference else None

        if not self.bam_path.exists():
            raise FileNotFoundError(f"BAM file not found: {self.bam_path}")

        # Check for index
        index_extensions = ['.bai', '.csi']
        has_index = any(
            self.bam_path.with_suffix(self.bam_path.suffix + ext).exists()
            or Path(str(self.bam_path) + ext).exists()
            for ext in index_extensions
        )
        if not has_index:
            raise FileNotFoundError(
                f"BAM index not found for: {self.bam_path}. "
                "Please index with 'samtools index'."
            )

        # Open file to check it's readable
        self._open_args = {"filename": str(self.bam_path)}
        if self.reference:
            self._open_args["reference_filename"] = str(self.reference)

    def _open_bam(self) -> pysam.AlignmentFile:
        """Open BAM file for reading."""
        return pysam.AlignmentFile(**self._open_args)

    @property
    def chromosomes(self) -> List[str]:
        """List of chromosome names in the BAM file."""
        with self._open_bam() as bam:
            return list(bam.references)

    @property
    def header(self) -> dict:
        """BAM header as dictionary."""
        with self._open_bam() as bam:
            return dict(bam.header)

    def count_reads(
        self,
        regions: List[GenomicRegion],
        count_method: str = "overlap"
    ) -> ReadCountResult:
        """Count reads in a list of genomic regions.

        Args:
            regions: List of genomic regions to count
            count_method: How to count reads:
                - "overlap": Count reads that overlap the region (default)
                - "start": Count reads whose start position is in the region
                - "midpoint": Count reads whose midpoint is in the region

        Returns:
            ReadCountResult with counts and statistics

        Raises:
            ValueError: If count_method is invalid
        """
        if not regions:
            return ReadCountResult(
                counts=np.array([], dtype=np.int32),
                regions=[],
                total_reads=0,
                filtered_reads=0
            )

        if count_method not in ("overlap", "start", "midpoint"):
            raise ValueError(f"Invalid count_method: {count_method}")

        counts = np.zeros(len(regions), dtype=np.int32)
        total_reads = 0
        filtered_reads = 0

        with self._open_bam() as bam:
            for i, region in enumerate(regions):
                count = 0

                try:
                    for read in bam.fetch(region.chrom, region.start, region.end):
                        total_reads += 1

                        # Filter by flags
                        if read.flag & self.EXCLUDE_FLAGS:
                            filtered_reads += 1
                            continue

                        # Filter by mapping quality
                        # Note: MAPQ 0 means multiple alignments, often included
                        if read.mapping_quality != 0 and read.mapping_quality < self.min_mapq:
                            filtered_reads += 1
                            continue

                        # Count based on method
                        if count_method == "overlap":
                            # pysam.fetch already filters by overlap
                            count += 1
                        elif count_method == "start":
                            if region.start <= read.reference_start < region.end:
                                count += 1
                        elif count_method == "midpoint":
                            midpoint = (read.reference_start + read.reference_end) // 2
                            if region.start <= midpoint < region.end:
                                count += 1
                except ValueError:
                    # Chromosome not in BAM
                    pass

                counts[i] = count

        return ReadCountResult(
            counts=counts,
            regions=regions,
            total_reads=total_reads,
            filtered_reads=filtered_reads
        )

    def count_reads_by_chromosome(
        self,
        regions: List[GenomicRegion],
        count_method: str = "overlap"
    ) -> Iterator[Tuple[str, ReadCountResult]]:
        """Count reads in regions, grouped by chromosome.

        This is more efficient than count_reads() when regions span
        multiple chromosomes, as it processes one chromosome at a time.

        Args:
            regions: List of genomic regions
            count_method: How to count reads (see count_reads)

        Yields:
            Tuples of (chromosome, ReadCountResult)
        """
        # Group regions by chromosome
        regions_by_chrom: dict = {}
        for region in regions:
            if region.chrom not in regions_by_chrom:
                regions_by_chrom[region.chrom] = []
            regions_by_chrom[region.chrom].append(region)

        # Process each chromosome
        for chrom in sorted(regions_by_chrom.keys()):
            chrom_regions = regions_by_chrom[chrom]
            result = self.count_reads(chrom_regions, count_method)
            yield chrom, result

    def count_reads_streaming(
        self,
        regions: List[GenomicRegion],
        chunk_size: int = 10000
    ) -> Iterator[Tuple[int, int, np.ndarray]]:
        """Count reads in chunks for memory efficiency.

        This method processes regions in chunks and yields partial results,
        suitable for very large region lists.

        Args:
            regions: List of genomic regions
            chunk_size: Number of regions per chunk

        Yields:
            Tuples of (start_index, end_index, counts_array)
        """
        for start in range(0, len(regions), chunk_size):
            end = min(start + chunk_size, len(regions))
            chunk_regions = regions[start:end]
            result = self.count_reads(chunk_regions)
            yield start, end, result.counts


def count_reads_in_regions(
    bam_path: Union[str, Path],
    regions: List[GenomicRegion],
    min_mapq: int = 20,
    reference: Optional[Union[str, Path]] = None,
    count_method: str = "overlap"
) -> ReadCountResult:
    """Convenience function for one-off read counting.

    Args:
        bam_path: Path to BAM/CRAM file
        regions: List of genomic regions
        min_mapq: Minimum mapping quality
        reference: Reference FASTA (for CRAM)
        count_method: How to count reads

    Returns:
        ReadCountResult with counts

    Example:
        >>> regions = [GenomicRegion("chr1", 1000, 2000, "gene1")]
        >>> result = count_reads_in_regions("sample.bam", regions)
        >>> print(result.counts)
    """
    reader = BAMReader(bam_path, min_mapq=min_mapq, reference=reference)
    return reader.count_reads(regions, count_method=count_method)


def load_regions_from_bed(
    bed_path: Union[str, Path],
    skip_header: bool = False
) -> List[GenomicRegion]:
    """Load genomic regions from a BED file.

    Args:
        bed_path: Path to BED file
        skip_header: Skip first line if it's a header

    Returns:
        List of GenomicRegion objects

    Note:
        BED format is 0-based, half-open [start, end)
    """
    regions = []
    bed_path = Path(bed_path)

    with open(bed_path) as f:
        for i, line in enumerate(f):
            if skip_header and i == 0:
                continue

            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track'):
                continue

            fields = line.split('\t')
            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else ""

            regions.append(GenomicRegion(chrom, start, end, name))

    return regions
