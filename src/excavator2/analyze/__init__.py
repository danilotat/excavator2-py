"""
EXCAVATOR2 analysis module.

This module provides high-level Python interfaces to the C++ HSLM
segmentation and FastCall CNV classification algorithms, as well as
the complete analysis pipeline.

Example:
    >>> from excavator2.analyze import CNVAnalyzer, AnalysisParameters
    >>> from excavator2.prepare import load_normalized_counts
    >>> # Load data
    >>> test = load_normalized_counts("test.NRC.h5")
    >>> control = load_normalized_counts("control.NRC.h5")
    >>> # Run analysis
    >>> params = AnalysisParameters(omega=0.1, cellularity=0.8)
    >>> analyzer = CNVAnalyzer(params)
    >>> result = analyzer.analyze_paired(test, control)
    >>> for cnv in result.get_cnvs():
    ...     print(f"{cnv.chrom}:{cnv.start}-{cnv.end} CN={cnv.absolute_cn}")
"""

# Segmentation
from excavator2.analyze.segment import (
    HSLMSegmenter,
    Segment,
    SegmentationResult,
    segment,
)

# CNV Classification
from excavator2.analyze.call import (
    FastCallCaller,
    CNVCall,
    ClassificationResult,
    CopyNumberState,
    call,
)

# Log2 ratio computation
from excavator2.analyze.ratio import (
    Log2RatioResult,
    compute_log2_ratio,
    compute_log2_ratio_pooled,
    apply_cellularity_correction,
)

# Analysis pipeline
from excavator2.analyze.pipeline import (
    CNVSegment,
    ChromosomeResult,
    AnalysisResult,
    AnalysisParameters,
    CNVAnalyzer,
    analyze_sample_pair,
)

__all__ = [
    # Segmentation classes and functions
    "HSLMSegmenter",
    "Segment",
    "SegmentationResult",
    "segment",
    # Classification classes and functions
    "FastCallCaller",
    "CNVCall",
    "ClassificationResult",
    "CopyNumberState",
    "call",
    # Log2 ratio computation
    "Log2RatioResult",
    "compute_log2_ratio",
    "compute_log2_ratio_pooled",
    "apply_cellularity_correction",
    # Analysis pipeline
    "CNVSegment",
    "ChromosomeResult",
    "AnalysisResult",
    "AnalysisParameters",
    "CNVAnalyzer",
    "analyze_sample_pair",
]
