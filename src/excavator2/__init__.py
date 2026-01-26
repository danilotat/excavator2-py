"""
EXCAVATOR2 - CNV detection from whole-exome sequencing data.

This is a complete rewrite of EXCAVATOR2 in Python/C++.
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
