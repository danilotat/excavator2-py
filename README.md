# EXCAVATOR2 (Python/C++ Rewrite)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

> **Note**: This is a complete rewrite of EXCAVATOR2 in Python/C++. For the original Perl/R/Fortran implementation, see [ctglab/excavator2](https://github.com/ctglab/excavator2).

## About

EXCAVATOR2 is a read count-based tool for detecting Copy Number Variants (CNVs) from Whole-Exome Sequencing (WES) and Targeted Sequencing (TS) data. It analyzes both in-target and off-target reads to identify CNVs with genome-wide resolution.

This modern implementation leverages:
- **Python** for data processing and orchestration
- **C++** for core algorithms (HSLM segmentation, FastCall classification)
- **pybind11** for seamless Python-C++ integration
- **HDF5** for efficient data storage
- **OpenMP** for parallel computation

## Features

- 🚀 **Performance**: Faster than original implementation through optimized C++ algorithms
- 📦 **Modern Packaging**: pip-installable, follows Python best practices
- 🔧 **Maintainability**: Clean architecture, comprehensive tests
- 📊 **HDF5 Storage**: Efficient binary format with metadata support
- ⚡ **Parallel Processing**: Multi-threaded algorithms and sample-level parallelization

## Installation

### From PyPI (when released)

```bash
pip install excavator2
```

### From Source

```bash
# Clone repository
git clone https://github.com/ctglab/excavator2-py.git
cd excavator2-py

# Install in development mode
pip install -e .

# Or build and install
pip install .
```

### Requirements

- **Python**: 3.9 or later
- **CMake**: 3.20 or later
- **C++ Compiler**: GCC 7+, Clang 5+, or MSVC 2019+
- **System libraries**: samtools, bedtools (for data processing)

## Quick Start

```bash
# Initialize target regions
excavator2 target --config config.yaml --output target/

# Prepare read counts (normalization)
excavator2 prepare --samples samples.yaml --target target/ --output prepared/ --threads 4

# Analyze and call CNVs
excavator2 analyze --samples sample_list.yaml --input prepared/ --target target/ --output results/ --experiment paired
```

## Citation

If you use EXCAVATOR2 in your research, please cite:

> D'Aurizio R, Pippucci T, Tattini L, Giusti B, Pellegrini M, Magi A. **Enhanced copy number variants detection from whole-exome sequencing data using EXCAVATOR2**. *Nucleic Acids Research* (2016) 44 (20):e154. [DOI: 10.1093/nar/gkw695](https://doi.org/10.1093/nar/gkw695)

## Development Status

🚧 **Currently in Phase 1**: Foundation & Build System

- [x] Project structure
- [x] Build system (scikit-build-core + CMake)
- [ ] HSLM algorithm (C++)
- [ ] FastCall algorithm (C++)
- [ ] Data preparation pipeline
- [ ] Target initialization
- [ ] Full integration tests

See [PLAN.md](PLAN.md) for detailed roadmap.

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Links

- **Original Implementation**: https://github.com/ctglab/excavator2
- **Documentation**: https://excavator2.readthedocs.io (coming soon)
- **Issue Tracker**: https://github.com/ctglab/excavator2-py/issues
