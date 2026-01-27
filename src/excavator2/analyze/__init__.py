"""
EXCAVATOR2 analysis module.

This module provides high-level Python interfaces to the C++ HSLM
segmentation and FastCall CNV classification algorithms.

Example:
    >>> from excavator2.analyze import HSLMSegmenter, FastCallCaller
    >>> # Segment log2 ratios
    >>> segmenter = HSLMSegmenter(omega=0.1)
    >>> seg_result = segmenter.segment(log2_ratios, positions)
    >>> # Classify segments
    >>> caller = FastCallCaller(cellularity=0.8)
    >>> call_result = caller.call([s.mean for s in seg_result.segments])
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
]
