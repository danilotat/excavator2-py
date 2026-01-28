"""
EXCAVATOR2 I/O module.

This module provides utilities for reading and writing various file formats
used in the EXCAVATOR2 pipeline.
"""

# BAM/CRAM reading
from excavator2.io.bam import (
    BAMReader,
    GenomicRegion,
    ReadCountResult,
    count_reads_in_regions,
    load_regions_from_bed,
)

# BED and coordinate file parsing
from excavator2.io.bed import (
    TargetRegion,
    ChromosomeInfo,
    GapRegion,
    load_target_bed,
    load_chromosome_coordinates,
    load_gap_file,
    regions_to_bed,
    merge_overlapping_regions,
    get_chromosome_regions,
    normalize_chromosome_name,
)

# BigWig reading (mappability)
from excavator2.io.bigwig import (
    BigWigReader,
    MappabilityResult,
    get_mappability_for_regions,
    check_bigwig_coverage,
)

# FASTA reading (reference genome)
from excavator2.io.fasta import (
    FastaReader,
    SequenceResult,
    GCResult,
    get_gc_content_for_regions,
    get_sequence_for_region,
)

# VCF/BED output is imported lazily to avoid circular imports
# Use: from excavator2.io.vcf import write_vcf

__all__ = [
    # BAM/CRAM
    "BAMReader",
    "GenomicRegion",
    "ReadCountResult",
    "count_reads_in_regions",
    "load_regions_from_bed",
    # BED
    "TargetRegion",
    "ChromosomeInfo",
    "GapRegion",
    "load_target_bed",
    "load_chromosome_coordinates",
    "load_gap_file",
    "regions_to_bed",
    "merge_overlapping_regions",
    "get_chromosome_regions",
    "normalize_chromosome_name",
    # BigWig
    "BigWigReader",
    "MappabilityResult",
    "get_mappability_for_regions",
    "check_bigwig_coverage",
    # FASTA
    "FastaReader",
    "SequenceResult",
    "GCResult",
    "get_gc_content_for_regions",
    "get_sequence_for_region",
    # VCF/BED output (lazy imports)
    "write_vcf",
    "write_vcf_regions",
    "write_vcf_windows",
    "write_bed",
    "write_segments_tsv",
    "write_fastcall_bed",
]


def __getattr__(name):
    """Lazy import for vcf module to avoid circular imports."""
    vcf_exports = {
        'write_vcf', 'write_vcf_regions', 'write_vcf_windows',
        'write_bed', 'write_segments_tsv', 'write_fastcall_bed'
    }
    if name in vcf_exports:
        from excavator2.io import vcf
        return getattr(vcf, name)
    raise AttributeError(f"module 'excavator2.io' has no attribute '{name}'")
