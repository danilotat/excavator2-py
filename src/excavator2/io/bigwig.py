"""
BigWig file reading utilities.

This module provides functions for reading BigWig files, particularly
for extracting mappability scores from mappability track files.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union, Tuple
import logging

import numpy as np

try:
    import pyBigWig

    HAS_PYBIGWIG = True
except ImportError:
    HAS_PYBIGWIG = False
    pyBigWig = None

logger = logging.getLogger(__name__)


@dataclass
class MappabilityResult:
    """Result of mappability extraction.

    Attributes:
        values: Mappability values per region (0.0-1.0)
        chromosomes: Chromosome names for each region
        starts: Start positions for each region
        ends: End positions for each region
        valid_mask: Boolean mask of regions with valid values
    """

    values: np.ndarray
    chromosomes: List[str]
    starts: np.ndarray
    ends: np.ndarray
    valid_mask: np.ndarray

    @property
    def n_regions(self) -> int:
        """Number of regions."""
        return len(self.values)

    @property
    def n_valid(self) -> int:
        """Number of regions with valid mappability values."""
        return int(np.sum(self.valid_mask))

    @property
    def mean_mappability(self) -> float:
        """Mean mappability of valid regions."""
        if self.n_valid == 0:
            return 0.0
        return float(np.mean(self.values[self.valid_mask]))


class BigWigReader:
    """BigWig file reader for mappability extraction.

    This class provides methods to extract mappability scores from BigWig
    files, typically used for k-mer mappability tracks.

    Example:
        >>> reader = BigWigReader("k100.Bismap.bw")
        >>> result = reader.get_mappability(["chr1", "chr1"], [1000, 2000], [1500, 2500])
        >>> print(f"Mean mappability: {result.mean_mappability:.3f}")

    Args:
        bigwig_path: Path to BigWig file
    """

    def __init__(self, bigwig_path: Union[str, Path]):
        if not HAS_PYBIGWIG:
            raise ImportError(
                "pyBigWig is required for BigWig reading. " "Install with: pip install pyBigWig"
            )

        self.bigwig_path = Path(bigwig_path)

        if not self.bigwig_path.exists():
            raise FileNotFoundError(f"BigWig file not found: {self.bigwig_path}")

        # Verify the file can be opened
        try:
            bw = pyBigWig.open(str(self.bigwig_path))
            self._chroms = bw.chroms()
            bw.close()
        except Exception as e:
            raise ValueError(f"Cannot open BigWig file {self.bigwig_path}: {e}")

    @property
    def chromosomes(self) -> List[str]:
        """List of chromosomes in the BigWig file."""
        return list(self._chroms.keys())

    @property
    def chromosome_sizes(self) -> dict:
        """Dictionary of chromosome sizes."""
        return dict(self._chroms)

    def get_mappability(
        self, chromosomes: List[str], starts: List[int], ends: List[int], stat: str = "mean"
    ) -> MappabilityResult:
        """Get mappability values for a list of regions.

        Args:
            chromosomes: List of chromosome names
            starts: List of start positions (0-based)
            ends: List of end positions (exclusive)
            stat: Statistic to compute over each region:
                - "mean": Mean mappability (default)
                - "min": Minimum mappability
                - "max": Maximum mappability

        Returns:
            MappabilityResult with values for each region
        """
        if len(chromosomes) != len(starts) or len(starts) != len(ends):
            raise ValueError("chromosomes, starts, and ends must have same length")

        n_regions = len(chromosomes)
        values = np.zeros(n_regions, dtype=np.float64)
        valid_mask = np.ones(n_regions, dtype=bool)

        bw = pyBigWig.open(str(self.bigwig_path))

        try:
            for i, (chrom, start, end) in enumerate(zip(chromosomes, starts, ends)):
                # Handle chromosome naming (with/without chr prefix)
                query_chrom = self._resolve_chromosome(chrom)

                if query_chrom is None:
                    logger.warning(f"Chromosome {chrom} not found in BigWig file")
                    values[i] = 0.0
                    valid_mask[i] = False
                    continue

                try:
                    # Get stats for the region
                    result = bw.stats(query_chrom, start, end, type=stat)
                    if result and result[0] is not None:
                        values[i] = result[0]
                    else:
                        values[i] = 0.0
                        valid_mask[i] = False
                except Exception as e:
                    logger.warning(f"Error getting mappability for {chrom}:{start}-{end}: {e}")
                    values[i] = 0.0
                    valid_mask[i] = False
        finally:
            bw.close()

        return MappabilityResult(
            values=values,
            chromosomes=chromosomes,
            starts=np.array(starts, dtype=np.int64),
            ends=np.array(ends, dtype=np.int64),
            valid_mask=valid_mask,
        )

    def get_values(self, chromosome: str, start: int, end: int) -> Optional[np.ndarray]:
        """Get base-level mappability values for a region.

        Args:
            chromosome: Chromosome name
            start: Start position (0-based)
            end: End position (exclusive)

        Returns:
            Array of mappability values (one per base), or None if error
        """
        query_chrom = self._resolve_chromosome(chromosome)
        if query_chrom is None:
            return None

        bw = pyBigWig.open(str(self.bigwig_path))
        try:
            values = bw.values(query_chrom, start, end)
            if values is None:
                return None
            return np.array(values, dtype=np.float64)
        except Exception as e:
            logger.warning(f"Error getting values for {chromosome}:{start}-{end}: {e}")
            return None
        finally:
            bw.close()

    def _resolve_chromosome(self, chrom: str) -> Optional[str]:
        """Resolve chromosome name to one that exists in the file.

        Args:
            chrom: Chromosome name to resolve

        Returns:
            Resolved chromosome name, or None if not found
        """
        if chrom in self._chroms:
            return chrom

        # Try alternate naming conventions
        if chrom.startswith("chr"):
            alt_chrom = chrom[3:]
        else:
            alt_chrom = f"chr{chrom}"

        if alt_chrom in self._chroms:
            return alt_chrom

        return None


def get_mappability_for_regions(
    bigwig_path: Union[str, Path],
    chromosomes: List[str],
    starts: List[int],
    ends: List[int],
    stat: str = "mean",
) -> np.ndarray:
    """Convenience function to get mappability values for regions.

    Args:
        bigwig_path: Path to BigWig file
        chromosomes: List of chromosome names
        starts: List of start positions (0-based)
        ends: List of end positions (exclusive)
        stat: Statistic to compute ("mean", "min", "max")

    Returns:
        Array of mappability values
    """
    reader = BigWigReader(bigwig_path)
    result = reader.get_mappability(chromosomes, starts, ends, stat=stat)
    return result.values


def check_bigwig_coverage(
    bigwig_path: Union[str, Path], chromosome: str, start: int, end: int, min_coverage: float = 0.5
) -> bool:
    """Check if a region has sufficient BigWig coverage.

    Args:
        bigwig_path: Path to BigWig file
        chromosome: Chromosome name
        start: Start position
        end: End position
        min_coverage: Minimum fraction of bases with values (0.0-1.0)

    Returns:
        True if region has sufficient coverage
    """
    reader = BigWigReader(bigwig_path)
    values = reader.get_values(chromosome, start, end)

    if values is None:
        return False

    # Count non-NaN values
    valid_count = np.sum(~np.isnan(values))
    coverage = valid_count / len(values)

    return coverage >= min_coverage
