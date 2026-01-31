"""
FASTA file reading utilities.

This module provides functions for reading reference genome FASTA files,
extracting sequences and computing GC content using pysam.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Optional, Union, Tuple
import logging

import numpy as np
import pysam

logger = logging.getLogger(__name__)


@dataclass
class SequenceResult:
    """Result of sequence extraction.

    Attributes:
        sequences: List of extracted sequences
        chromosomes: Chromosome names for each sequence
        starts: Start positions for each sequence
        ends: End positions for each sequence
    """

    sequences: List[str]
    chromosomes: List[str]
    starts: np.ndarray
    ends: np.ndarray

    @property
    def n_sequences(self) -> int:
        """Number of sequences."""
        return len(self.sequences)


@dataclass
class GCResult:
    """Result of GC content calculation.

    Attributes:
        gc_content: GC fraction for each region (0.0-1.0)
        chromosomes: Chromosome names for each region
        starts: Start positions for each region
        ends: End positions for each region
        n_bases: Number of bases in each region (excluding N)
    """

    gc_content: np.ndarray
    chromosomes: List[str]
    starts: np.ndarray
    ends: np.ndarray
    n_bases: np.ndarray

    @property
    def n_regions(self) -> int:
        """Number of regions."""
        return len(self.gc_content)

    @property
    def mean_gc(self) -> float:
        """Mean GC content."""
        if len(self.gc_content) == 0:
            return 0.0
        return float(np.mean(self.gc_content))


class FastaReader:
    """FASTA file reader for sequence extraction and GC calculation.

    This class provides methods to extract sequences and compute GC content
    from reference genome FASTA files using pysam.

    Example:
        >>> reader = FastaReader("GRCh38.fasta")
        >>> gc = reader.get_gc_content(["chr1"], [1000], [2000])
        >>> print(f"GC content: {gc.gc_content[0]:.2%}")

    Args:
        fasta_path: Path to FASTA file (must be indexed with samtools faidx)
    """

    def __init__(self, fasta_path: Union[str, Path]):
        self.fasta_path = Path(fasta_path)

        if not self.fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {self.fasta_path}")

        # Check for FASTA index
        fai_path = Path(str(self.fasta_path) + ".fai")
        if not fai_path.exists():
            logger.warning(
                f"FASTA index not found: {fai_path}. "
                "Creating index with pysam (may take time for large files)..."
            )
            pysam.faidx(str(self.fasta_path))

        # Open and get chromosome list
        try:
            self._fasta = pysam.FastaFile(str(self.fasta_path))
            self._chroms = dict(zip(self._fasta.references, self._fasta.lengths))
        except Exception as e:
            raise ValueError(f"Cannot open FASTA file {self.fasta_path}: {e}")

    def close(self) -> None:
        """Close the FASTA file."""
        if hasattr(self, "_fasta") and self._fasta:
            self._fasta.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    @property
    def chromosomes(self) -> List[str]:
        """List of chromosome names in the FASTA file."""
        return list(self._chroms.keys())

    @property
    def chromosome_sizes(self) -> Dict[str, int]:
        """Dictionary of chromosome sizes."""
        return dict(self._chroms)

    def get_sequence(self, chromosome: str, start: int, end: int) -> Optional[str]:
        """Get the sequence for a genomic region.

        Args:
            chromosome: Chromosome name
            start: Start position (0-based)
            end: End position (exclusive)

        Returns:
            Sequence string (uppercase), or None if error
        """
        query_chrom = self._resolve_chromosome(chromosome)
        if query_chrom is None:
            logger.warning(f"Chromosome {chromosome} not found in FASTA file")
            return None

        try:
            # pysam uses 0-based coordinates
            seq = self._fasta.fetch(query_chrom, start, end)
            return seq.upper()
        except Exception as e:
            logger.warning(f"Error fetching {chromosome}:{start}-{end}: {e}")
            return None

    def get_sequences(
        self, chromosomes: List[str], starts: List[int], ends: List[int]
    ) -> SequenceResult:
        """Get sequences for multiple regions.

        Args:
            chromosomes: List of chromosome names
            starts: List of start positions (0-based)
            ends: List of end positions (exclusive)

        Returns:
            SequenceResult with sequences for each region
        """
        if len(chromosomes) != len(starts) or len(starts) != len(ends):
            raise ValueError("chromosomes, starts, and ends must have same length")

        sequences = []
        for chrom, start, end in zip(chromosomes, starts, ends):
            seq = self.get_sequence(chrom, start, end)
            sequences.append(seq if seq else "")

        return SequenceResult(
            sequences=sequences,
            chromosomes=chromosomes,
            starts=np.array(starts, dtype=np.int64),
            ends=np.array(ends, dtype=np.int64),
        )

    def get_gc_content(
        self, chromosomes: List[str], starts: List[int], ends: List[int]
    ) -> GCResult:
        """Calculate GC content for multiple regions.

        Args:
            chromosomes: List of chromosome names
            starts: List of start positions (0-based)
            ends: List of end positions (exclusive)

        Returns:
            GCResult with GC fraction for each region (0.0-1.0)
        """
        if len(chromosomes) != len(starts) or len(starts) != len(ends):
            raise ValueError("chromosomes, starts, and ends must have same length")

        n_regions = len(chromosomes)
        gc_content = np.zeros(n_regions, dtype=np.float64)
        n_bases = np.zeros(n_regions, dtype=np.int64)

        for i, (chrom, start, end) in enumerate(zip(chromosomes, starts, ends)):
            seq = self.get_sequence(chrom, start, end)
            if seq:
                gc, total = self._compute_gc(seq)
                gc_content[i] = gc
                n_bases[i] = total
            else:
                gc_content[i] = 0.0
                n_bases[i] = 0

        return GCResult(
            gc_content=gc_content,
            chromosomes=chromosomes,
            starts=np.array(starts, dtype=np.int64),
            ends=np.array(ends, dtype=np.int64),
            n_bases=n_bases,
        )

    def get_base_at_position(self, chromosome: str, position: int) -> Optional[str]:
        """Get the reference base at a specific position.

        Args:
            chromosome: Chromosome name
            position: Position (0-based)

        Returns:
            Single-character base (uppercase), or None if error
        """
        query_chrom = self._resolve_chromosome(chromosome)
        if query_chrom is None:
            return None

        try:
            # Fetch single base
            base = self._fasta.fetch(query_chrom, position, position + 1)
            return base.upper()
        except Exception:
            return None

    def get_bases_at_positions(self, chromosomes: List[str], positions: List[int]) -> List[str]:
        """Get reference bases at multiple positions.

        Args:
            chromosomes: List of chromosome names
            positions: List of positions (0-based)

        Returns:
            List of single-character bases
        """
        if len(chromosomes) != len(positions):
            raise ValueError("chromosomes and positions must have same length")

        bases = []
        for chrom, pos in zip(chromosomes, positions):
            base = self.get_base_at_position(chrom, pos)
            bases.append(base if base else "N")

        return bases

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

    @staticmethod
    def _compute_gc(sequence: str) -> Tuple[float, int]:
        """Compute GC fraction for a sequence.

        Args:
            sequence: DNA sequence (uppercase)

        Returns:
            Tuple of (gc_fraction, total_valid_bases)
            GC fraction is 0.0 if no valid bases
        """
        # Count bases (excluding N and other ambiguous bases)
        g_count = sequence.count("G")
        c_count = sequence.count("C")
        a_count = sequence.count("A")
        t_count = sequence.count("T")

        total = g_count + c_count + a_count + t_count

        if total == 0:
            return 0.0, 0

        gc_fraction = (g_count + c_count) / total
        return gc_fraction, total


def get_gc_content_for_regions(
    fasta_path: Union[str, Path], chromosomes: List[str], starts: List[int], ends: List[int]
) -> np.ndarray:
    """Convenience function to get GC content for regions.

    Args:
        fasta_path: Path to FASTA file
        chromosomes: List of chromosome names
        starts: List of start positions (0-based)
        ends: List of end positions (exclusive)

    Returns:
        Array of GC fractions (0.0-1.0)
    """
    with FastaReader(fasta_path) as reader:
        result = reader.get_gc_content(chromosomes, starts, ends)
    return result.gc_content


def get_sequence_for_region(
    fasta_path: Union[str, Path], chromosome: str, start: int, end: int
) -> Optional[str]:
    """Convenience function to get sequence for a single region.

    Args:
        fasta_path: Path to FASTA file
        chromosome: Chromosome name
        start: Start position (0-based)
        end: End position (exclusive)

    Returns:
        Sequence string, or None if error
    """
    with FastaReader(fasta_path) as reader:
        return reader.get_sequence(chromosome, start, end)
