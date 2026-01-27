"""
EXCAVATOR2 data preparation module.

This module provides functionality for:
- Reading and counting reads from BAM/CRAM files
- Normalizing read counts for systematic biases (size, mappability, GC)
- Saving/loading processed data in HDF5 format

Example:
    >>> from excavator2.prepare import ReadCountProcessor, ReadCountNormalizer
    >>> # Count reads
    >>> processor = ReadCountProcessor("target.h5")
    >>> sample_data = processor.process_sample("sample.bam", "sample1")
    >>> # Normalize
    >>> normalizer = ReadCountNormalizer()
    >>> result = normalizer.normalize(sample_data)
"""

# Read counting
from excavator2.prepare.readcount import (
    WindowData,
    SampleReadCounts,
    ReadCountProcessor,
    save_read_counts,
    load_read_counts,
)

# Normalization
from excavator2.prepare.normalize import (
    NormalizationResult,
    ReadCountNormalizer,
    normalize_read_counts,
    save_normalized_counts,
    load_normalized_counts,
)

__all__ = [
    # Read counting
    "WindowData",
    "SampleReadCounts",
    "ReadCountProcessor",
    "save_read_counts",
    "load_read_counts",
    # Normalization
    "NormalizationResult",
    "ReadCountNormalizer",
    "normalize_read_counts",
    "save_normalized_counts",
    "load_normalized_counts",
]
