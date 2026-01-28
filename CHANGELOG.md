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

#### Added - Phase 4: Python Algorithm Wrappers
- High-level Python wrappers for HSLM and FastCall algorithms
  - `src/excavator2/analyze/segment.py` - HSLM segmentation wrapper
  - `src/excavator2/analyze/call.py` - FastCall classification wrapper
  - `src/excavator2/analyze/__init__.py` - Module exports
- HSLMSegmenter class:
  - Pythonic interface with parameter validation
  - `segment()` method for single-chromosome segmentation
  - `segment_multi()` method for multi-sample joint segmentation
  - Returns `SegmentationResult` with `Segment` objects containing genomic coordinates
  - Input validation (empty arrays, NaN/Inf, length mismatches)
- FastCallCaller class:
  - Pythonic interface with parameter validation
  - `call()` method for CNV classification
  - Returns `ClassificationResult` with `CNVCall` objects
  - `CopyNumberState` enum with labels and helper methods
  - Helper methods: `get_cnvs()`, `get_deletions()`, `get_gains()`
- Convenience functions:
  - `segment()` for one-off HSLM segmentation
  - `call()` for one-off FastCall classification
- Data classes with full type hints:
  - `Segment`: start/end indices, genomic positions, mean, n_probes
  - `SegmentationResult`: segments, breakpoints, state_path, success/error
  - `CNVCall`: cn_call, absolute_cn, state, probability, segment_mean
  - `ClassificationResult`: calls, state_means/sds/priors, EM stats
- Integration tests for full HSLM→FastCall pipeline

#### Notes - Phase 4
- All tests passing (92/92)
- 40 new tests for Python wrappers
- 73% code coverage for analyze module (91% for call.py, 95% for segment.py)
- End-to-end pipeline integration verified

#### Added - Phase 5: Data Preparation Pipeline
- BAM/CRAM file I/O module (`src/excavator2/io/bam.py`)
  - `BAMReader` class for reading BAM/CRAM files with pysam
  - `GenomicRegion` dataclass for genomic intervals
  - `ReadCountResult` dataclass for read counting results
  - Read counting with MAPQ filtering and duplicate exclusion
  - `load_regions_from_bed()` for BED file parsing
  - Support for CRAM files with reference FASTA
- Read count processing module (`src/excavator2/prepare/readcount.py`)
  - `WindowData` dataclass for window metadata (GC, mappability, class)
  - `SampleReadCounts` dataclass for raw read counts per sample
  - `ReadCountProcessor` class for BAM→counts pipeline
  - Chromosome-wise data access methods
  - IN-target vs OFF-target region masks
  - HDF5 save/load functions for read counts
- Normalization module (`src/excavator2/prepare/normalize.py`)
  - `ReadCountNormalizer` class implementing median-based corrections
  - `NormalizationResult` dataclass for normalized counts
  - Three-stage normalization pipeline:
    1. Size/length normalization (5 bp bins)
    2. Mappability normalization (5% bins)
    3. GC-content normalization (5% bins)
  - Separate handling for IN-target and OFF-target regions
  - Zero replacement with minimum non-zero value
  - HDF5 save/load functions for normalized counts
- CLI `prepare` command fully implemented
  - Sample sheet YAML parsing
  - Read counting from BAM files
  - Normalization with progress output
  - Output to HDF5 format (`.RC.h5` and `.NRC.h5`)
  - Support for parallel processing (via `--threads`)
  - Reference FASTA option for CRAM files

#### Notes - Phase 5
- All tests passing (122/122)
- 30 new tests for prepare module
- 67% overall code coverage
- Normalization algorithm matches original R implementation
- HDF5 format for efficient storage and cross-platform compatibility

#### Added - Phase 6: Target Initialization
- BED file parsing module (`src/excavator2/io/bed.py`)
  - `TargetRegion` dataclass for target regions
  - `ChromosomeInfo` dataclass for chromosome metadata
  - `GapRegion` dataclass for gap annotations
  - `load_target_bed()` for BED file parsing with filtering
  - `load_chromosome_coordinates()` for chromosome size files
  - `load_gap_file()` for UCSC gap file parsing
  - `merge_overlapping_regions()` for region merging
  - Chromosome name normalization utilities
- BigWig reading module (`src/excavator2/io/bigwig.py`)
  - `BigWigReader` class for mappability extraction
  - `MappabilityResult` dataclass for results
  - Support for pyBigWig library
  - Mean/min/max statistics per region
- FASTA reading module (`src/excavator2/io/fasta.py`)
  - `FastaReader` class for reference genome access
  - `GCResult` dataclass for GC content results
  - `SequenceResult` dataclass for sequence extraction
  - Auto-indexing with pysam if .fai missing
- Target filtering module (`src/excavator2/target/filter.py`)
  - `AnalysisWindow` dataclass for analysis windows
  - `create_analysis_windows()` for IN-target and OUT-target windows
  - `filter_gap_overlapping_windows()` for gap filtering
  - `create_filtered_target()` main entry point
  - Window renumbering utilities
- GC content module (`src/excavator2/target/gc_content.py`)
  - `calculate_gc_content()` for window GC calculation
  - `calculate_gc_content_by_chromosome()` for per-chromosome processing
  - `get_window_gc_stats()` for statistics
- Mappability module (`src/excavator2/target/mappability.py`)
  - `extract_mappability()` for window mappability extraction
  - `extract_mappability_by_chromosome()` for per-chromosome processing
  - `get_mappability_stats()` for statistics
  - `filter_low_mappability_windows()` for quality filtering
- Target initialization module (`src/excavator2/target/init.py`)
  - `TargetData` dataclass for complete target data
  - `initialize_target()` main initialization function
  - `save_target_data()` HDF5 serialization
  - `load_target_data()` HDF5 deserialization
- CLI `target` command fully implemented
  - Config YAML parsing (Reference + Target sections)
  - Window creation from BED regions
  - Gap/centromere filtering
  - GC content calculation from FASTA
  - Mappability extraction from BigWig
  - Output to HDF5 format
  - Settings file generation

#### Notes - Phase 6
- All tests passing (161/161)
- 39 new tests for target module
- 61% overall code coverage
- Window creation matches original R/Perl implementation
- HDF5 format compatible with prepare module

#### Added - Phase 7: Integration & Full Pipeline
- Log2 ratio computation module (`src/excavator2/analyze/ratio.py`)
  - `Log2RatioResult` dataclass with chromosome/region accessors
  - `compute_log2_ratio()` for paired test vs control
  - `compute_log2_ratio_pooled()` for pooled control mode
  - `apply_cellularity_correction()` for tumor purity adjustment
  - Median centering with separate IN/OUT target handling
- Analysis pipeline module (`src/excavator2/analyze/pipeline.py`)
  - `CNVSegment` dataclass with genomic coordinates and CN calls
  - `ChromosomeResult` for per-chromosome results
  - `AnalysisResult` for complete analysis output
  - `AnalysisParameters` for HSLM + FastCall configuration
  - `CNVAnalyzer` class orchestrating full pipeline:
    - `analyze_paired()` for matched test/control samples
    - `analyze_pooled()` for pooled control experiments
  - Segment filtering by minimum probe count
- VCF/BED output module (`src/excavator2/io/vcf.py`)
  - `write_vcf()` - VCF 4.2 format with SV annotations
  - `write_vcf_regions()` - CNVs only for downstream analysis
  - `write_vcf_windows()` - All segments for visualization
  - `write_bed()` - BED format with CNV type and probability
  - `write_segments_tsv()` - Full segment details (HSLM format)
  - `write_fastcall_bed()` - FastCall format with CN calls
- CLI `analyze` command fully implemented
  - Sample file list YAML parsing (T1/C1, T2/C2 format)
  - Paired mode: matched test/control pairs
  - Pooled mode: test samples vs pooled controls
  - Parameters YAML with HSLM and FastCall settings
  - Multi-sample batch processing
  - Per-sample output directories
  - Output formats: VCF, BED, TSV, FastCall BED
  - Analysis settings saved with each run

#### Notes - Phase 7
- All tests passing (192/192)
- 31 new tests for Phase 7 integration
- 62% overall code coverage
- Full end-to-end pipeline integration verified
- Both paired and pooled experimental modes working
- Output formats match original EXCAVATOR2 conventions

## Roadmap

- ~~**Phase 2 (Weeks 3-5)**: C++ HSLM Algorithm~~ COMPLETED
- ~~**Phase 3 (Weeks 6-8)**: C++ FastCall Algorithm~~ COMPLETED
- ~~**Phase 4 (Week 9)**: Python Algorithm Wrappers~~ COMPLETED
- ~~**Phase 5 (Weeks 10-12)**: Data Preparation Pipeline~~ COMPLETED
- ~~**Phase 6 (Weeks 13-14)**: Target Initialization~~ COMPLETED
- ~~**Phase 7 (Weeks 15-16)**: Integration & Full Pipeline~~ COMPLETED
- **Phase 8 (Week 17)**: Documentation, CI/CD, Release
