"""
Read count processing and WMRC calculation.

This module handles reading and processing raw read counts from BAM files,
organizing them by chromosome, and preparing them for normalization.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union
import logging

import numpy as np
import h5py

from excavator2.io.bam import BAMReader, GenomicRegion, ReadCountResult

logger = logging.getLogger(__name__)


def _count_chromosome_reads(args):
    """Count reads for a single chromosome (picklable for multiprocessing).

    Args:
        args: Tuple of (chrom, bam_path, min_mapq, reference, regions)
            where regions is a list of (start, end, name) tuples

    Returns:
        Tuple of (counts_array, total_reads, filtered_reads)
    """
    chrom, bam_path, min_mapq, reference, regions = args

    from excavator2.io.bam import BAMReader, GenomicRegion

    # Create fresh BAM reader for this process
    reader = BAMReader(bam_path, min_mapq=min_mapq, reference=reference)

    # Convert to GenomicRegion objects
    genomic_regions = [
        GenomicRegion(chrom=chrom, start=start, end=end, name=name) for start, end, name in regions
    ]

    # Count reads
    result = reader.count_reads(genomic_regions, count_method="overlap")

    return result.counts, result.total_reads, result.filtered_reads


@dataclass
class WindowData:
    """Data for a single genomic window/exon.

    Attributes:
        chrom: Chromosome name
        start: Start position (0-based)
        end: End position (exclusive)
        position: Midpoint position
        gc_content: GC content (0-1)
        mappability: Mappability score (0-1)
        region_class: 'IN' for in-target, 'OUT' for off-target
        name: Optional region name (e.g., gene name)
    """

    chrom: str
    start: int
    end: int
    position: int
    gc_content: float
    mappability: float
    region_class: str  # 'IN' or 'OUT'
    name: str = ""

    @property
    def length(self) -> int:
        """Window length in base pairs."""
        return self.end - self.start


@dataclass
class SampleReadCounts:
    """Read counts for a single sample.

    Attributes:
        sample_name: Sample identifier
        raw_counts: Raw read counts per window
        windows: Window metadata
        chromosomes: List of unique chromosomes
    """

    sample_name: str
    raw_counts: np.ndarray
    windows: List[WindowData]
    chromosomes: List[str] = field(default_factory=list)

    def __post_init__(self):
        if not self.chromosomes:
            self.chromosomes = sorted(set(w.chrom for w in self.windows))

    @property
    def n_windows(self) -> int:
        """Total number of windows."""
        return len(self.windows)

    @property
    def total_reads(self) -> int:
        """Total read count across all windows."""
        return int(np.sum(self.raw_counts))

    def get_chromosome_mask(self, chrom: str) -> np.ndarray:
        """Get boolean mask for windows on a chromosome."""
        return np.array([w.chrom == chrom for w in self.windows])

    def get_chromosome_data(self, chrom: str) -> tuple:
        """Get counts and windows for a specific chromosome.

        Returns:
            Tuple of (counts, windows) for the chromosome
        """
        mask = self.get_chromosome_mask(chrom)
        counts = self.raw_counts[mask]
        windows = [w for w, m in zip(self.windows, mask) if m]
        return counts, windows

    def get_in_target_mask(self) -> np.ndarray:
        """Get boolean mask for in-target windows."""
        return np.array([w.region_class == "IN" for w in self.windows])

    def get_off_target_mask(self) -> np.ndarray:
        """Get boolean mask for off-target windows."""
        return np.array([w.region_class == "OUT" for w in self.windows])


class ReadCountProcessor:
    """Process BAM files to generate read counts.

    This class handles the complete workflow of counting reads from BAM
    files in genomic windows defined by a target file.

    Example:
        >>> processor = ReadCountProcessor(target_path="target.h5")
        >>> sample_data = processor.process_sample("sample.bam", "sample1")
        >>> print(f"Total reads: {sample_data.total_reads}")

    Args:
        target_path: Path to HDF5 target file with window definitions
        min_mapq: Minimum mapping quality for read filtering
        reference: Path to reference FASTA (for CRAM files)
    """

    def __init__(
        self,
        target_path: Union[str, Path],
        min_mapq: int = 20,
        reference: Optional[Union[str, Path]] = None,
    ):
        self.target_path = Path(target_path)
        self.min_mapq = min_mapq
        self.reference = reference

        if not self.target_path.exists():
            raise FileNotFoundError(f"Target file not found: {self.target_path}")

        # Load window definitions from target file
        self.windows = self._load_target_windows()
        self.chromosomes = sorted(set(w.chrom for w in self.windows))
        logger.info(f"Loaded {len(self.windows)} windows from {len(self.chromosomes)} chromosomes")

    def _load_target_windows(self) -> List[WindowData]:
        """Load window definitions from HDF5 target file."""
        windows = []

        with h5py.File(self.target_path, "r") as f:
            # Check if we have the expected structure
            if "chromosomes" not in f:
                raise ValueError(f"Invalid target file format: missing 'chromosomes' group")

            for chrom in f["chromosomes"]:
                chrom_grp = f["chromosomes"][chrom]

                # Load arrays
                starts = chrom_grp["start"][:]
                ends = chrom_grp["end"][:]
                gc_content = chrom_grp["gc_content"][:]
                mappability = chrom_grp["mappability"][:]

                # Handle region class - may be bytes or string
                region_class = chrom_grp["class"][:]
                if region_class.dtype.kind == "S":  # byte string
                    region_class = [c.decode("utf-8") for c in region_class]
                else:
                    region_class = list(region_class)

                # Optional: gene names
                if "name" in chrom_grp:
                    names = chrom_grp["name"][:]
                    if names.dtype.kind == "S":
                        names = [n.decode("utf-8") for n in names]
                else:
                    names = [""] * len(starts)

                # Create WindowData objects
                for i in range(len(starts)):
                    windows.append(
                        WindowData(
                            chrom=chrom,
                            start=int(starts[i]),
                            end=int(ends[i]),
                            position=(int(starts[i]) + int(ends[i])) // 2,
                            gc_content=float(gc_content[i]),
                            mappability=float(mappability[i]),
                            region_class=region_class[i],
                            name=names[i] if i < len(names) else "",
                        )
                    )

        return windows

    def _windows_to_regions(self) -> List[GenomicRegion]:
        """Convert WindowData to GenomicRegion for BAM counting."""
        return [
            GenomicRegion(chrom=w.chrom, start=w.start, end=w.end, name=w.name)
            for w in self.windows
        ]

    def process_sample(
        self, bam_path: Union[str, Path], sample_name: str, n_threads: int = 1
    ) -> SampleReadCounts:
        """Process a BAM file and count reads in all windows.

        Args:
            bam_path: Path to BAM/CRAM file
            sample_name: Sample identifier
            n_threads: Number of parallel threads for chromosome processing

        Returns:
            SampleReadCounts with raw counts and metadata
        """
        logger.info(f"Processing sample: {sample_name}")
        logger.info(f"BAM file: {bam_path}")

        if n_threads > 1 and len(self.chromosomes) > 1:
            return self._process_sample_parallel(bam_path, sample_name, n_threads)

        # Sequential processing
        # Create BAM reader
        reader = BAMReader(bam_path, min_mapq=self.min_mapq, reference=self.reference)

        # Convert windows to regions
        regions = self._windows_to_regions()

        # Count reads
        logger.info(f"Counting reads in {len(regions)} windows...")
        result = reader.count_reads(regions, count_method="overlap")

        logger.info(f"Total reads processed: {result.total_reads}")
        logger.info(f"Filtered reads: {result.filtered_reads}")
        logger.info(f"Mean count per window: {result.mean_count:.2f}")

        return SampleReadCounts(
            sample_name=sample_name,
            raw_counts=result.counts,
            windows=self.windows.copy(),
            chromosomes=self.chromosomes.copy(),
        )

    def _process_sample_parallel(
        self, bam_path: Union[str, Path], sample_name: str, n_threads: int
    ) -> SampleReadCounts:
        """Process a BAM file with chromosome-level parallelization.

        Args:
            bam_path: Path to BAM/CRAM file
            sample_name: Sample identifier
            n_threads: Number of parallel threads

        Returns:
            SampleReadCounts with raw counts and metadata
        """
        from concurrent.futures import ProcessPoolExecutor, as_completed

        # Group windows by chromosome
        windows_by_chrom = {}
        for i, w in enumerate(self.windows):
            if w.chrom not in windows_by_chrom:
                windows_by_chrom[w.chrom] = []
            windows_by_chrom[w.chrom].append((i, w))

        # Prepare work items for each chromosome
        work_items = []
        for chrom in self.chromosomes:
            if chrom in windows_by_chrom:
                chrom_windows = windows_by_chrom[chrom]
                work_items.append(
                    (
                        chrom,
                        str(bam_path),
                        self.min_mapq,
                        str(self.reference) if self.reference else None,
                        [(w.start, w.end, w.name) for _, w in chrom_windows],
                    )
                )

        # Process chromosomes in parallel
        logger.info(f"Processing {len(work_items)} chromosomes with {n_threads} threads...")

        chrom_results = {}
        with ProcessPoolExecutor(max_workers=n_threads) as executor:
            futures = {
                executor.submit(_count_chromosome_reads, item): item[0] for item in work_items
            }

            for future in as_completed(futures):
                chrom = futures[future]
                try:
                    chrom_results[chrom] = future.result()
                    logger.debug(f"Completed chromosome {chrom}")
                except Exception as e:
                    logger.error(f"Failed to process chromosome {chrom}: {e}")
                    raise

        # Reassemble results in original window order
        counts = np.zeros(len(self.windows), dtype=np.int32)
        total_reads = 0
        filtered_reads = 0

        for chrom, (chrom_counts, chrom_total, chrom_filtered) in chrom_results.items():
            chrom_windows = windows_by_chrom[chrom]
            for (orig_idx, _), count in zip(chrom_windows, chrom_counts):
                counts[orig_idx] = count
            total_reads += chrom_total
            filtered_reads += chrom_filtered

        logger.info(f"Total reads processed: {total_reads}")
        logger.info(f"Filtered reads: {filtered_reads}")
        logger.info(f"Mean count per window: {counts.mean():.2f}")

        return SampleReadCounts(
            sample_name=sample_name,
            raw_counts=counts,
            windows=self.windows.copy(),
            chromosomes=self.chromosomes.copy(),
        )

    def process_samples_parallel(
        self, samples: Dict[str, Union[str, Path]], n_workers: int = 1
    ) -> Dict[str, SampleReadCounts]:
        """Process multiple samples, optionally in parallel.

        Args:
            samples: Dictionary mapping sample names to BAM paths
            n_workers: Number of parallel workers

        Returns:
            Dictionary mapping sample names to SampleReadCounts
        """
        from concurrent.futures import ProcessPoolExecutor, as_completed

        results = {}

        if n_workers <= 1:
            # Sequential processing
            for sample_name, bam_path in samples.items():
                results[sample_name] = self.process_sample(bam_path, sample_name)
        else:
            # Parallel processing
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                futures = {
                    executor.submit(self.process_sample, bam_path, sample_name): sample_name
                    for sample_name, bam_path in samples.items()
                }

                for future in as_completed(futures):
                    sample_name = futures[future]
                    try:
                        results[sample_name] = future.result()
                    except Exception as e:
                        logger.error(f"Failed to process {sample_name}: {e}")
                        raise

        return results


def save_read_counts(
    sample_data: SampleReadCounts, output_path: Union[str, Path], compression: str = "gzip"
) -> None:
    """Save read counts to HDF5 file.

    Args:
        sample_data: SampleReadCounts to save
        output_path: Path to output HDF5 file
        compression: Compression algorithm (default: gzip)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(output_path, "w") as f:
        # Metadata
        f.attrs["sample_name"] = sample_data.sample_name
        f.attrs["n_windows"] = sample_data.n_windows
        f.attrs["total_reads"] = sample_data.total_reads

        # Create chromosomes group
        chrom_grp = f.create_group("chromosomes")

        for chrom in sample_data.chromosomes:
            counts, windows = sample_data.get_chromosome_data(chrom)

            chr_grp = chrom_grp.create_group(chrom)

            # Store arrays
            chr_grp.create_dataset("raw_counts", data=counts, compression=compression)
            chr_grp.create_dataset(
                "start", data=[w.start for w in windows], compression=compression
            )
            chr_grp.create_dataset("end", data=[w.end for w in windows], compression=compression)
            chr_grp.create_dataset(
                "position", data=[w.position for w in windows], compression=compression
            )
            chr_grp.create_dataset(
                "gc_content", data=[w.gc_content for w in windows], compression=compression
            )
            chr_grp.create_dataset(
                "mappability", data=[w.mappability for w in windows], compression=compression
            )

            # String arrays need special handling
            region_class = np.array([w.region_class for w in windows], dtype="S10")
            chr_grp.create_dataset("class", data=region_class, compression=compression)

            names = np.array([w.name for w in windows], dtype=h5py.special_dtype(vlen=str))
            chr_grp.create_dataset("name", data=names, compression=compression)


def load_read_counts(input_path: Union[str, Path]) -> SampleReadCounts:
    """Load read counts from HDF5 file.

    Args:
        input_path: Path to HDF5 file

    Returns:
        SampleReadCounts object
    """
    input_path = Path(input_path)

    with h5py.File(input_path, "r") as f:
        sample_name = f.attrs["sample_name"]

        windows = []
        all_counts = []
        chromosomes = []

        for chrom in f["chromosomes"]:
            chromosomes.append(chrom)
            chr_grp = f["chromosomes"][chrom]

            counts = chr_grp["raw_counts"][:]
            starts = chr_grp["start"][:]
            ends = chr_grp["end"][:]
            positions = chr_grp["position"][:]
            gc_content = chr_grp["gc_content"][:]
            mappability = chr_grp["mappability"][:]

            region_class = chr_grp["class"][:]
            if region_class.dtype.kind == "S":
                region_class = [c.decode("utf-8") for c in region_class]

            names = chr_grp["name"][:]
            if hasattr(names, "astype"):
                names = [str(n) for n in names]

            all_counts.extend(counts)

            for i in range(len(counts)):
                windows.append(
                    WindowData(
                        chrom=chrom,
                        start=int(starts[i]),
                        end=int(ends[i]),
                        position=int(positions[i]),
                        gc_content=float(gc_content[i]),
                        mappability=float(mappability[i]),
                        region_class=region_class[i],
                        name=names[i] if i < len(names) else "",
                    )
                )

        return SampleReadCounts(
            sample_name=sample_name,
            raw_counts=np.array(all_counts, dtype=np.int32),
            windows=windows,
            chromosomes=sorted(chromosomes),
        )
