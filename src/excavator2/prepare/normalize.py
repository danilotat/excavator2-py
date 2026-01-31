"""
Read count normalization algorithms.

This module implements the normalization pipeline for correcting systematic
biases in read counts:
1. Size/length normalization
2. Mappability normalization
3. GC-content normalization

All corrections use median-based scaling within feature bins.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union
import logging

import numpy as np
import h5py

from excavator2.prepare.readcount import SampleReadCounts, WindowData

logger = logging.getLogger(__name__)


@dataclass
class NormalizationResult:
    """Result of read count normalization.

    Attributes:
        sample_name: Sample identifier
        normalized_counts: Normalized read counts
        windows: Window metadata
        chromosomes: List of unique chromosomes
        normalization_stats: Statistics from normalization steps
    """

    sample_name: str
    normalized_counts: np.ndarray
    windows: List[WindowData]
    chromosomes: List[str]
    normalization_stats: dict

    @property
    def n_windows(self) -> int:
        """Total number of windows."""
        return len(self.windows)


def _median_normalize(
    counts: np.ndarray, feature: np.ndarray, bin_size: float, min_bin_count: int = 1
) -> Tuple[np.ndarray, dict]:
    """Apply median-based normalization within feature bins.

    This is the core normalization algorithm used for size, mappability,
    and GC-content corrections. For each bin of the feature:
    - Calculate the median count in that bin
    - Scale all counts in the bin by (master_median / bin_median)

    Bin boundaries follow the original R implementation:
    - First bin: [bin_start, bin_end] (inclusive both ends)
    - Other bins: (bin_start, bin_end] (exclusive left, inclusive right)

    Args:
        counts: Read counts to normalize
        feature: Feature values to bin by (same length as counts)
        bin_size: Size of each feature bin
        min_bin_count: Minimum windows in a bin to apply correction (default: 1)

    Returns:
        Tuple of (normalized_counts, stats_dict)
    """
    if len(counts) == 0:
        return counts.copy(), {}

    # Work with float copy
    normalized = counts.astype(np.float64)

    # Calculate master median INCLUDING zeros (matching R's na.rm=T behavior)
    # R: MasterMedian <- median(RC, na.rm = T)
    # This includes zeros but excludes NaN values
    master_median = np.nanmedian(normalized)

    # If master median is 0, no normalization is possible
    if master_median == 0:
        return normalized, {"master_median": 0.0, "n_bins": 0}

    # Create bins
    feature_min = np.floor(np.min(feature) / bin_size) * bin_size
    feature_max = np.ceil(np.max(feature) / bin_size) * bin_size
    bins = np.arange(feature_min, feature_max + bin_size, bin_size)

    n_bins_corrected = 0
    bin_medians = []

    for i in range(len(bins) - 1):
        bin_start = bins[i]
        bin_end = bins[i + 1]

        # Find windows in this bin
        # Match R boundary handling:
        # - First bin: >= bin_start AND <= bin_end
        # - Other bins: > bin_start AND <= bin_end
        if i == 0:
            bin_mask = (feature >= bin_start) & (feature <= bin_end)
        else:
            bin_mask = (feature > bin_start) & (feature <= bin_end)

        bin_counts = normalized[bin_mask]

        # R only requires length(ind) > 0, so min_bin_count defaults to 1
        if len(bin_counts) < min_bin_count:
            continue

        # Calculate bin median INCLUDING zeros (matching R's na.rm=T)
        # R: m <- median(RC[ind], na.rm = T)
        bin_median = np.nanmedian(bin_counts)
        bin_medians.append(bin_median)

        # R checks: if (m > 0) - skip normalization if bin median is 0
        # This is critical: R leaves data unchanged for bins with median=0
        if bin_median > 0:
            # Scale counts in this bin
            scale_factor = master_median / bin_median
            normalized[bin_mask] = normalized[bin_mask] * scale_factor
            n_bins_corrected += 1

    stats = {
        "master_median": float(master_median),
        "n_bins": n_bins_corrected,
        "bin_size": bin_size,
        "bin_medians": bin_medians,
    }

    return normalized, stats


class ReadCountNormalizer:
    """Normalize read counts for systematic biases.

    This class implements the full normalization pipeline:
    1. Size/length normalization (corrects for exon length)
    2. Mappability normalization (corrects for mappability differences)
    3. GC-content normalization (corrects for GC bias)

    In-target and off-target regions are normalized separately,
    as they have different characteristics.

    Example:
        >>> normalizer = ReadCountNormalizer()
        >>> result = normalizer.normalize(sample_data)
        >>> print(f"Normalized {result.n_windows} windows")

    Args:
        size_bin: Bin size for length normalization (bp). Default: 5
        mappability_bin: Bin size for mappability (%). Default: 5
        gc_bin: Bin size for GC content (%). Default: 5
        min_bin_count: Minimum windows per bin. Default: 1 (matches original R)
    """

    def __init__(
        self,
        size_bin: float = 5.0,
        mappability_bin: float = 5.0,
        gc_bin: float = 5.0,
        min_bin_count: int = 1,
    ):
        self.size_bin = size_bin
        self.mappability_bin = mappability_bin
        self.gc_bin = gc_bin
        self.min_bin_count = min_bin_count

    def normalize(
        self,
        sample_data: SampleReadCounts,
        normalize_size: bool = True,
        normalize_mappability: bool = True,
        normalize_gc: bool = True,
    ) -> NormalizationResult:
        """Apply full normalization pipeline.

        Args:
            sample_data: Raw read counts to normalize
            normalize_size: Apply size normalization (default: True)
            normalize_mappability: Apply mappability normalization (default: True)
            normalize_gc: Apply GC normalization (default: True)

        Returns:
            NormalizationResult with normalized counts
        """
        logger.info(f"Normalizing sample: {sample_data.sample_name}")

        # Extract arrays
        raw_counts = sample_data.raw_counts.astype(np.float64)
        lengths = np.array([w.length for w in sample_data.windows])

        # CRITICAL: Convert raw counts to WMRC (Width-normalized Mapped Read Count)
        # by dividing by window length. This is essential before any other normalization.
        # Original R: RCTMatrixL <- t(t(RCTL) / L)
        # Avoid division by zero for zero-length windows
        safe_lengths = np.maximum(lengths, 1.0)
        counts = raw_counts / safe_lengths
        logger.info(f"  Converted to WMRC (counts/length)")
        mappability = np.array([w.mappability * 100 for w in sample_data.windows])  # Convert to %
        gc_content = np.array([w.gc_content * 100 for w in sample_data.windows])  # Convert to %

        # Get masks for IN and OUT target regions
        in_mask = sample_data.get_in_target_mask()
        out_mask = sample_data.get_off_target_mask()

        logger.info(f"  IN-target windows: {np.sum(in_mask)}")
        logger.info(f"  OFF-target windows: {np.sum(out_mask)}")

        # Initialize normalized counts
        normalized = counts.copy()
        stats = {}

        # Process IN-target regions
        if np.any(in_mask):
            normalized[in_mask], in_stats = self._normalize_subset(
                counts[in_mask],
                lengths[in_mask],
                mappability[in_mask],
                gc_content[in_mask],
                normalize_size=normalize_size,
                normalize_mappability=normalize_mappability,
                normalize_gc=normalize_gc,
                label="IN",
            )
            stats["in_target"] = in_stats

        # Process OFF-target regions (skip size normalization)
        if np.any(out_mask):
            normalized[out_mask], out_stats = self._normalize_subset(
                counts[out_mask],
                lengths[out_mask],
                mappability[out_mask],
                gc_content[out_mask],
                normalize_size=False,  # OFF-target doesn't use size normalization
                normalize_mappability=normalize_mappability,
                normalize_gc=normalize_gc,
                label="OUT",
            )
            stats["off_target"] = out_stats

        # Replace zeros with minimum non-zero value
        nonzero = normalized[normalized > 0]
        if len(nonzero) > 0:
            min_nonzero = np.min(nonzero)
            normalized[normalized == 0] = min_nonzero
            stats["min_value_replacement"] = float(min_nonzero)

        return NormalizationResult(
            sample_name=sample_data.sample_name,
            normalized_counts=normalized,
            windows=sample_data.windows,
            chromosomes=sample_data.chromosomes,
            normalization_stats=stats,
        )

    def _normalize_subset(
        self,
        counts: np.ndarray,
        lengths: np.ndarray,
        mappability: np.ndarray,
        gc_content: np.ndarray,
        normalize_size: bool,
        normalize_mappability: bool,
        normalize_gc: bool,
        label: str,
    ) -> Tuple[np.ndarray, dict]:
        """Normalize a subset of windows (IN or OUT target).

        Args:
            counts: Read counts
            lengths: Window lengths
            mappability: Mappability scores (0-100%)
            gc_content: GC content (0-100%)
            normalize_size: Apply size normalization
            normalize_mappability: Apply mappability normalization
            normalize_gc: Apply GC normalization
            label: Label for logging

        Returns:
            Tuple of (normalized_counts, stats)
        """
        normalized = counts.astype(np.float64)
        stats = {}

        # Step 1: Size/length normalization
        if normalize_size:
            logger.info(f"  {label}: Applying size normalization (bin={self.size_bin}bp)")
            normalized, size_stats = _median_normalize(
                normalized, lengths, self.size_bin, self.min_bin_count
            )
            stats["size"] = size_stats
            logger.info(
                f"    Corrected {size_stats['n_bins']} bins, master median={size_stats['master_median']:.2f}"
            )

        # Step 2: Mappability normalization
        if normalize_mappability:
            logger.info(
                f"  {label}: Applying mappability normalization (bin={self.mappability_bin}%)"
            )
            normalized, map_stats = _median_normalize(
                normalized, mappability, self.mappability_bin, self.min_bin_count
            )
            stats["mappability"] = map_stats
            logger.info(
                f"    Corrected {map_stats['n_bins']} bins, master median={map_stats['master_median']:.2f}"
            )

        # Step 3: GC-content normalization
        if normalize_gc:
            logger.info(f"  {label}: Applying GC normalization (bin={self.gc_bin}%)")
            normalized, gc_stats = _median_normalize(
                normalized, gc_content, self.gc_bin, self.min_bin_count
            )
            stats["gc"] = gc_stats
            logger.info(
                f"    Corrected {gc_stats['n_bins']} bins, master median={gc_stats['master_median']:.2f}"
            )

        return normalized, stats


def normalize_read_counts(
    sample_data: SampleReadCounts,
    size_bin: float = 5.0,
    mappability_bin: float = 5.0,
    gc_bin: float = 5.0,
) -> NormalizationResult:
    """Convenience function for normalizing read counts.

    Args:
        sample_data: Raw read counts
        size_bin: Bin size for length normalization (bp)
        mappability_bin: Bin size for mappability (%)
        gc_bin: Bin size for GC content (%)

    Returns:
        NormalizationResult with normalized counts

    Example:
        >>> result = normalize_read_counts(sample_data)
        >>> print(f"Normalized mean: {np.mean(result.normalized_counts):.2f}")
    """
    normalizer = ReadCountNormalizer(
        size_bin=size_bin, mappability_bin=mappability_bin, gc_bin=gc_bin
    )
    return normalizer.normalize(sample_data)


def save_normalized_counts(
    result: NormalizationResult, output_path: Union[str, Path], compression: str = "gzip"
) -> None:
    """Save normalized counts to HDF5 file.

    Args:
        result: NormalizationResult to save
        output_path: Path to output HDF5 file
        compression: Compression algorithm
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(output_path, "w") as f:
        # Metadata
        f.attrs["sample_name"] = result.sample_name
        f.attrs["n_windows"] = result.n_windows

        # Create chromosomes group
        chrom_grp = f.create_group("chromosomes")

        for chrom in result.chromosomes:
            # Get data for this chromosome
            mask = np.array([w.chrom == chrom for w in result.windows])
            counts = result.normalized_counts[mask]
            windows = [w for w, m in zip(result.windows, mask) if m]

            chr_grp = chrom_grp.create_group(chrom)

            # Store arrays
            chr_grp.create_dataset("normalized_counts", data=counts, compression=compression)
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

            # String arrays
            region_class = np.array([w.region_class for w in windows], dtype="S10")
            chr_grp.create_dataset("class", data=region_class, compression=compression)

            names = np.array([w.name for w in windows], dtype=h5py.special_dtype(vlen=str))
            chr_grp.create_dataset("name", data=names, compression=compression)

        # Store normalization stats
        stats_grp = f.create_group("normalization_stats")
        _save_nested_dict(stats_grp, result.normalization_stats)


def _save_nested_dict(group: h5py.Group, data: dict) -> None:
    """Recursively save nested dictionary to HDF5 group."""
    for key, value in data.items():
        if isinstance(value, dict):
            subgroup = group.create_group(str(key))
            _save_nested_dict(subgroup, value)
        elif isinstance(value, (list, np.ndarray)):
            try:
                group.create_dataset(str(key), data=np.array(value))
            except (TypeError, ValueError):
                # Skip non-numeric lists
                pass
        elif isinstance(value, (int, float)):
            group.attrs[str(key)] = value
        elif isinstance(value, str):
            group.attrs[str(key)] = value


def load_normalized_counts(input_path: Union[str, Path]) -> NormalizationResult:
    """Load normalized counts from HDF5 file.

    Args:
        input_path: Path to HDF5 file

    Returns:
        NormalizationResult object
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

            counts = chr_grp["normalized_counts"][:]
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

        # Load normalization stats (simplified)
        stats = {}
        if "normalization_stats" in f:
            stats = _load_nested_dict(f["normalization_stats"])

        return NormalizationResult(
            sample_name=sample_name,
            normalized_counts=np.array(all_counts),
            windows=windows,
            chromosomes=sorted(chromosomes),
            normalization_stats=stats,
        )


def _load_nested_dict(group: h5py.Group) -> dict:
    """Recursively load nested dictionary from HDF5 group."""
    result = {}

    # Load attributes
    for key, value in group.attrs.items():
        result[key] = value

    # Load datasets and subgroups
    for key in group.keys():
        item = group[key]
        if isinstance(item, h5py.Group):
            result[key] = _load_nested_dict(item)
        elif isinstance(item, h5py.Dataset):
            result[key] = item[:]

    return result
