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

## Roadmap

- **Phase 2 (Weeks 3-5)**: C++ HSLM Algorithm
- **Phase 3 (Weeks 6-8)**: C++ FastCall Algorithm
- **Phase 4 (Week 9)**: Python Algorithm Wrappers
- **Phase 5 (Weeks 10-12)**: Data Preparation Pipeline
- **Phase 6 (Weeks 13-14)**: Target Initialization
- **Phase 7 (Weeks 15-16)**: Integration & Full Pipeline
- **Phase 8 (Week 17)**: Documentation, CI/CD, Release
