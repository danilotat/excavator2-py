# Changelog

All notable changes to the EXCAVATOR2 Python/C++ rewrite will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### [3.0.0-alpha] - 2026-01-26

#### Added - Phase 1: Foundation & Build System
- Initial project structure with src/ layout (PEP 420)
- Python package structure (`src/excavator2/`)
- C++ source structure (`cpp/include`, `cpp/src`)
- Build system configuration:
  - `pyproject.toml` with scikit-build-core + pybind11
  - `CMakeLists.txt` for C++ compilation
  - OpenMP support detection and linking
- Minimal pybind11 "hello world" C++ module (`_excavator_core`)
- CLI skeleton with click:
  - `excavator2 target` (stub)
  - `excavator2 prepare` (stub)
  - `excavator2 analyze` (stub)
- Testing infrastructure:
  - pytest configuration
  - Basic unit tests
  - Test coverage reporting
- Development files:
  - README.md with project overview
  - .gitignore
  - CHANGELOG.md

#### Technical Details
- Python version: 3.9+
- C++ standard: C++17
- Build system: scikit-build-core 0.11+ with CMake 3.20+
- Bindings: pybind11 3.0+
- OpenMP: Enabled (detected at build time)

#### Notes
- All tests passing (4/4)
- C++ module successfully compiles and loads
- CLI commands are stubs awaiting implementation in future phases

#### Added - Phase 2: C++ HSLM Algorithm
- HSLM (Heterogeneous Shifting Level Model) segmentation algorithm in C++
  - `cpp/include/excavator/common.hpp` - Log-space utilities
  - `cpp/include/excavator/hslm.hpp` - HSLM algorithm interface
  - `cpp/src/common.cpp` - Utility implementations
  - `cpp/src/hslm.cpp` - Full HSLM implementation
- Core algorithm components ported from Fortran:
  - `elnsum()` - Numerically stable log-sum-exp
  - `logsumexp()` - Log-sum-exp for vectors
  - `TRANSEMISI` - HMM transition/emission matrix computation
  - `BIOVITERBII` - Viterbi algorithm for optimal state path
- Python bindings via pybind11:
  - `excavator2._excavator_core.hslm` module
  - `HSLMParameters` class for algorithm configuration
  - `HSLMResult` class for segmentation results
  - `HSLM` class with `segment()` and `segment_multi()` methods
  - Convenience `segment()` function for single-call usage
- Features:
  - Distance-dependent transition probabilities
  - Multi-sample segmentation support
  - Configurable parameters (omega, theta, step_eta, n_states)
  - Segment filtering by minimum size
  - OpenMP parallelization for emission matrix computation

#### Notes - Phase 2
- All tests passing (27/27)
- HSLM correctly detects CNV segments in synthetic data
- Algorithm matches expected behavior from original Fortran/R implementation

#### Added - Phase 3: C++ FastCall Algorithm
- FastCall CNV classification algorithm in C++
  - `cpp/include/excavator/fastcall.hpp` - FastCall algorithm interface
  - `cpp/src/fastcall.cpp` - Full FastCall implementation
- Core algorithm components ported from R:
  - `truncated_gaussian_pdf()` - Truncated Gaussian PDF (gfct)
  - `run_em()` - Expectation-Maximization algorithm
  - `e_step()` - E-step using truncated Gaussians
  - `m_step()` - M-step updating parameters
  - `compute_posteriors()` - Posterior probability calculation
  - `assign_labels()` - CN state assignment from posteriors
- Python bindings via pybind11:
  - `excavator2._excavator_core.fastcall` module
  - `FastCallParameters` class for algorithm configuration
  - `SegmentCall` class for individual segment results
  - `FastCallResult` class for classification results
  - `FastCall` class with `call()` method
  - Convenience `call()` function for single-call usage
- Features:
  - 5-state Gaussian mixture model (CN: 0, 1, 2, 3, 4+)
  - EM algorithm with convergence detection
  - Configurable parameters (cellularity, thrd, thru)
  - Posterior probability for each call
  - Handles edge cases (outliers, constant data, extreme values)

#### Notes - Phase 3
- All tests passing (52/52)
- FastCall correctly classifies CN states in synthetic data
- EM converges in 1-10 iterations for typical data
- Algorithm matches expected behavior from original R implementation

## Roadmap

- ~~**Phase 2 (Weeks 3-5)**: C++ HSLM Algorithm~~ COMPLETED
- ~~**Phase 3 (Weeks 6-8)**: C++ FastCall Algorithm~~ COMPLETED
- **Phase 4 (Week 9)**: Python Algorithm Wrappers
- **Phase 5 (Weeks 10-12)**: Data Preparation Pipeline
- **Phase 6 (Weeks 13-14)**: Target Initialization
- **Phase 7 (Weeks 15-16)**: Integration & Full Pipeline
- **Phase 8 (Week 17)**: Documentation, CI/CD, Release
