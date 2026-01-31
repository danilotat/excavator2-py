"""
Target initialization module.

This module provides the main target initialization functionality,
combining window creation, gap filtering, and feature extraction
(GC content, mappability) into a single pipeline.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Optional, Union
import logging
from datetime import datetime

import numpy as np
import h5py

from excavator2.target.filter import (
    AnalysisWindow,
    FilteredTargetResult,
    create_filtered_target,
    renumber_windows,
)
from excavator2.target.gc_content import calculate_gc_content, get_window_gc_stats
from excavator2.target.mappability import extract_mappability, get_mappability_stats

logger = logging.getLogger(__name__)


@dataclass
class TargetData:
    """Complete target data for CNV analysis.

    Attributes:
        windows: List of analysis windows
        gc_content: GC fraction for each window (0.0-1.0)
        mappability: Mappability score for each window (0.0-1.0)
        chromosomes: List of unique chromosomes
        target_name: Name of the target
        assembly: Reference assembly name (e.g., hg38)
        window_size: Window size used
        metadata: Additional metadata
    """

    windows: List[AnalysisWindow]
    gc_content: np.ndarray
    mappability: np.ndarray
    chromosomes: List[str]
    target_name: str
    assembly: str
    window_size: int
    metadata: dict = field(default_factory=dict)

    @property
    def n_windows(self) -> int:
        """Total number of windows."""
        return len(self.windows)

    @property
    def n_in_target(self) -> int:
        """Number of IN-target windows."""
        return sum(1 for w in self.windows if w.region_class == "IN")

    @property
    def n_out_target(self) -> int:
        """Number of OUT-target windows."""
        return sum(1 for w in self.windows if w.region_class == "OUT")

    def get_chromosome_data(self, chromosome: str) -> dict:
        """Get data for a specific chromosome.

        Args:
            chromosome: Chromosome name

        Returns:
            Dictionary with windows, gc_content, mappability for the chromosome
        """
        # Build mask for this chromosome
        chrom_variants = {chromosome}
        if chromosome.startswith("chr"):
            chrom_variants.add(chromosome[3:])
        else:
            chrom_variants.add(f"chr{chromosome}")

        mask = np.array([w.chrom in chrom_variants for w in self.windows])
        indices = np.where(mask)[0]

        return {
            "windows": [self.windows[i] for i in indices],
            "gc_content": self.gc_content[mask],
            "mappability": self.mappability[mask],
            "n_windows": int(np.sum(mask)),
        }


def initialize_target(
    bed_path: Union[str, Path],
    fasta_path: Union[str, Path],
    bigwig_path: Union[str, Path],
    chromosome_path: Union[str, Path],
    gap_path: Union[str, Path],
    target_name: str,
    assembly: str,
    window_size: int,
    flank: int = 200,
    chromosomes: Optional[List[str]] = None,
) -> TargetData:
    """Initialize target data for CNV analysis.

    This is the main entry point for target preparation. It:
    1. Creates analysis windows from target BED regions
    2. Filters windows overlapping genomic gaps
    3. Calculates GC content from reference FASTA
    4. Extracts mappability from BigWig file

    Args:
        bed_path: Path to target BED file
        fasta_path: Path to reference genome FASTA
        bigwig_path: Path to mappability BigWig file
        chromosome_path: Path to chromosome coordinates file
        gap_path: Path to gap annotations file
        target_name: Name for this target set
        assembly: Reference assembly name (e.g., hg38)
        window_size: Size of analysis windows
        flank: Buffer around gaps (default: 200bp)
        chromosomes: If provided, only process these chromosomes

    Returns:
        TargetData with windows and features
    """
    logger.info(f"Initializing target: {target_name}")
    logger.info(f"Assembly: {assembly}, Window size: {window_size}")

    # Step 1: Create filtered windows
    logger.info("Step 1: Creating analysis windows...")
    filtered_result = create_filtered_target(
        bed_path=bed_path,
        chromosome_path=chromosome_path,
        gap_path=gap_path,
        window_size=window_size,
        target_name=target_name,
        flank=flank,
        chromosomes=chromosomes,
    )

    # Renumber windows sequentially
    windows = renumber_windows(filtered_result.windows)

    # Step 2: Calculate GC content
    logger.info("Step 2: Calculating GC content...")
    gc_content = calculate_gc_content(windows, fasta_path)
    gc_stats = get_window_gc_stats(gc_content, windows)
    logger.info(
        f"  Mean GC: {gc_stats['mean']:.2%}, Range: {gc_stats['min']:.2%}-{gc_stats['max']:.2%}"
    )

    # Step 3: Extract mappability
    logger.info("Step 3: Extracting mappability...")
    mappability = extract_mappability(windows, bigwig_path)
    map_stats = get_mappability_stats(mappability, windows)
    logger.info(
        f"  Mean MAP: {map_stats['mean']:.3f}, Low MAP regions: {map_stats['low_mappability_count']}"
    )

    # Build metadata
    metadata = {
        "creation_date": datetime.now().isoformat(),
        "bed_path": str(bed_path),
        "fasta_path": str(fasta_path),
        "bigwig_path": str(bigwig_path),
        "flank": flank,
        "n_filtered_gaps": filtered_result.n_filtered,
        "gc_stats": gc_stats,
        "mappability_stats": map_stats,
    }

    logger.info(f"Target initialization complete: {len(windows)} windows")

    return TargetData(
        windows=windows,
        gc_content=gc_content,
        mappability=mappability,
        chromosomes=filtered_result.chromosomes,
        target_name=target_name,
        assembly=assembly,
        window_size=window_size,
        metadata=metadata,
    )


def save_target_data(
    target: TargetData, output_path: Union[str, Path], compression: str = "gzip"
) -> None:
    """Save target data to HDF5 file.

    Args:
        target: TargetData to save
        output_path: Output HDF5 file path
        compression: Compression algorithm (default: gzip)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Saving target data to {output_path}")

    with h5py.File(output_path, "w") as f:
        # Global metadata
        f.attrs["target_name"] = target.target_name
        f.attrs["assembly"] = target.assembly
        f.attrs["window_size"] = target.window_size
        f.attrs["n_windows"] = target.n_windows
        f.attrs["n_in_target"] = target.n_in_target
        f.attrs["n_out_target"] = target.n_out_target

        # Store metadata dict
        if target.metadata:
            meta_grp = f.create_group("metadata")
            _save_metadata(meta_grp, target.metadata)

        # Store by chromosome
        chrom_grp = f.create_group("chromosomes")

        for chrom in target.chromosomes:
            chrom_data = target.get_chromosome_data(chrom)
            chr_grp = chrom_grp.create_group(chrom)

            windows = chrom_data["windows"]
            n = len(windows)

            # Window coordinates
            chr_grp.create_dataset(
                "start", data=[w.start for w in windows], compression=compression
            )
            chr_grp.create_dataset("end", data=[w.end for w in windows], compression=compression)
            chr_grp.create_dataset(
                "position", data=[w.position for w in windows], compression=compression
            )

            # Window metadata
            window_ids = np.array(
                [w.window_id for w in windows], dtype=h5py.special_dtype(vlen=str)
            )
            chr_grp.create_dataset("window_id", data=window_ids, compression=compression)

            region_class = np.array([w.region_class for w in windows], dtype="S3")
            chr_grp.create_dataset("class", data=region_class, compression=compression)

            # Features
            chr_grp.create_dataset(
                "gc_content", data=chrom_data["gc_content"], compression=compression
            )
            chr_grp.create_dataset(
                "mappability", data=chrom_data["mappability"], compression=compression
            )

            chr_grp.attrs["n_windows"] = n

    logger.info(
        f"Saved target data: {target.n_windows} windows across {len(target.chromosomes)} chromosomes"
    )


def _save_metadata(group: h5py.Group, metadata: dict) -> None:
    """Save metadata dictionary to HDF5 group."""
    for key, value in metadata.items():
        if isinstance(value, dict):
            # Create subgroup for nested dicts
            subgroup = group.create_group(str(key))
            _save_metadata(subgroup, value)
        elif isinstance(value, (int, float)):
            group.attrs[str(key)] = value
        elif isinstance(value, str):
            group.attrs[str(key)] = value
        elif isinstance(value, (list, np.ndarray)):
            try:
                group.create_dataset(str(key), data=np.array(value))
            except (TypeError, ValueError):
                # Skip non-numeric lists
                pass


def load_target_data(input_path: Union[str, Path]) -> TargetData:
    """Load target data from HDF5 file.

    Args:
        input_path: Path to HDF5 file

    Returns:
        TargetData object
    """
    input_path = Path(input_path)

    with h5py.File(input_path, "r") as f:
        target_name = f.attrs["target_name"]
        assembly = f.attrs["assembly"]
        window_size = int(f.attrs["window_size"])

        # Load metadata
        metadata = {}
        if "metadata" in f:
            metadata = _load_metadata(f["metadata"])

        # Load windows and features
        windows = []
        gc_list = []
        map_list = []
        chromosomes = []

        for chrom in f["chromosomes"]:
            chromosomes.append(chrom)
            chr_grp = f["chromosomes"][chrom]

            starts = chr_grp["start"][:]
            ends = chr_grp["end"][:]
            positions = chr_grp["position"][:]
            window_ids = chr_grp["window_id"][:]
            region_class = chr_grp["class"][:]
            gc_content = chr_grp["gc_content"][:]
            mappability = chr_grp["mappability"][:]

            # Decode bytes
            if region_class.dtype.kind == "S":
                region_class = [c.decode("utf-8") for c in region_class]
            if hasattr(window_ids, "astype"):
                window_ids = [str(w) for w in window_ids]

            for i in range(len(starts)):
                windows.append(
                    AnalysisWindow(
                        chrom=chrom,
                        start=int(starts[i]),
                        end=int(ends[i]),
                        window_id=window_ids[i] if i < len(window_ids) else f"a{i}",
                        region_class=region_class[i],
                        position=int(positions[i]),
                    )
                )
                gc_list.append(gc_content[i])
                map_list.append(mappability[i])

        return TargetData(
            windows=windows,
            gc_content=np.array(gc_list),
            mappability=np.array(map_list),
            chromosomes=sorted(chromosomes),
            target_name=target_name,
            assembly=assembly,
            window_size=window_size,
            metadata=metadata,
        )


def _load_metadata(group: h5py.Group) -> dict:
    """Load metadata dictionary from HDF5 group."""
    result = {}

    # Load attributes
    for key, value in group.attrs.items():
        result[key] = value

    # Load datasets and subgroups
    for key in group.keys():
        item = group[key]
        if isinstance(item, h5py.Group):
            result[key] = _load_metadata(item)
        elif isinstance(item, h5py.Dataset):
            result[key] = item[:]

    return result
