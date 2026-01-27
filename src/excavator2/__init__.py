"""
EXCAVATOR2 - CNV detection from whole-exome sequencing data.

This is a complete rewrite of EXCAVATOR2 in Python/C++.

Example:
    >>> from excavator2.analyze import HSLMSegmenter, FastCallCaller
    >>> segmenter = HSLMSegmenter()
    >>> seg_result = segmenter.segment(log2_ratios, positions)
    >>> caller = FastCallCaller()
    >>> call_result = caller.call([s.mean for s in seg_result.segments])
"""

__version__ = "3.0.0"
__author__ = "EXCAVATOR2 Development Team"

# Make sure this gets set when the package is imported
try:
    from excavator2._excavator_core import __version__ as _cpp_version
    if _cpp_version != __version__:
        import warnings
        warnings.warn(f"Version mismatch: Python={__version__}, C++={_cpp_version}")
except ImportError:
    # C++ module not yet built
    pass

# Expose analyze module for convenience
from excavator2 import analyze

__all__ = [
    "__version__",
    "__author__",
    "analyze",
]
