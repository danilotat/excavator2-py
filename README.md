# EXCAVATOR2 (Python/C++ Rewrite)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Build](https://github.com/ctglab/excavator2-py/actions/workflows/build.yml/badge.svg)](https://github.com/ctglab/excavator2-py/actions/workflows/build.yml)
[![Test](https://github.com/ctglab/excavator2-py/actions/workflows/test.yml/badge.svg)](https://github.com/ctglab/excavator2-py/actions/workflows/test.yml)
[![Format](https://github.com/ctglab/excavator2-py/actions/workflows/format.yml/badge.svg)](https://github.com/ctglab/excavator2-py/actions/workflows/format.yml)
[![Lint](https://github.com/ctglab/excavator2-py/actions/workflows/lint.yml/badge.svg)](https://github.com/ctglab/excavator2-py/actions/workflows/lint.yml)

## About

EXCAVATOR2 is a read count-based tool for detecting Copy Number Variants (CNVs) from Whole-Exome Sequencing (WES) and Targeted Sequencing (TS) data. It analyzes both in-target and off-target reads to identify CNVs with genome-wide resolution.

This is a complete rewrite of EXCAVATOR2 in Python/C++. For the original Perl/R/Fortran implementation, see [ctglab/excavator2](https://github.com/ctglab/excavator2).

### Implementation Status

This rewrite is fully functional and produces CNV calls. However, it is not yet 100% numerically identical to the original implementation. Minor differences may exist in edge cases due to floating-point handling and algorithmic refinements during porting. Validation against the original is ongoing.

The implementation uses:
- **Python** for data processing and orchestration
- **C++** for core algorithms (HSLM segmentation, FastCall classification)
- **pybind11** for Python-C++ integration
- **HDF5** for data storage

## Installation

### From Source

```bash
git clone https://github.com/ctglab/excavator2-py.git
cd excavator2-py
pip install -e .
```

### Requirements

- Python 3.9+
- CMake 3.20+
- C++ compiler with C++17 support (GCC 7+, Clang 5+, MSVC 2019+)

## Usage

```bash
# Initialize target regions
excavator2 target --config config.yaml --output target/

# Prepare read counts
excavator2 prepare --samples samples.yaml --target target/ --output prepared/ --threads 4

# Call CNVs
excavator2 analyze --samples sample_list.yaml --input prepared/ --target target/ --output results/ --experiment paired
```

## Citation

If you use EXCAVATOR2 in your research, please cite:

> D'Aurizio R, Pippucci T, Tattini L, Giusti B, Pellegrini M, Magi A. **Enhanced copy number variants detection from whole-exome sequencing data using EXCAVATOR2**. *Nucleic Acids Research* (2016) 44 (20):e154. [DOI: 10.1093/nar/gkw695](https://doi.org/10.1093/nar/gkw695)

## License

GNU General Public License v3.0. See [LICENSE](LICENSE).

## Links

- Original implementation: https://github.com/ctglab/excavator2
- Issue tracker: https://github.com/ctglab/excavator2-py/issues
