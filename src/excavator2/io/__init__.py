"""
EXCAVATOR2 I/O module.

This module provides utilities for reading and writing various file formats
used in the EXCAVATOR2 pipeline.
"""

from excavator2.io.bam import (
    BAMReader,
    GenomicRegion,
    ReadCountResult,
    count_reads_in_regions,
    load_regions_from_bed,
)

__all__ = [
    "BAMReader",
    "GenomicRegion",
    "ReadCountResult",
    "count_reads_in_regions",
    "load_regions_from_bed",
]
