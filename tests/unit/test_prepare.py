"""
Unit tests for EXCAVATOR2 data preparation module.

Tests the BAM I/O, read counting, and normalization functionality.
"""

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import h5py

# =============================================================================
# Tests for io/bam.py - GenomicRegion and ReadCountResult
# =============================================================================


class TestGenomicRegion:
    """Tests for GenomicRegion dataclass."""

    def test_creation(self):
        """Test basic GenomicRegion creation."""
        from excavator2.io.bam import GenomicRegion

        region = GenomicRegion("chr1", 1000, 2000, "gene1")

        assert region.chrom == "chr1"
        assert region.start == 1000
        assert region.end == 2000
        assert region.name == "gene1"

    def test_length_property(self):
        """Test length calculation."""
        from excavator2.io.bam import GenomicRegion

        region = GenomicRegion("chr1", 1000, 2000)
        assert region.length == 1000

        region2 = GenomicRegion("chr1", 0, 500)
        assert region2.length == 500

    def test_midpoint_property(self):
        """Test midpoint calculation."""
        from excavator2.io.bam import GenomicRegion

        region = GenomicRegion("chr1", 1000, 2000)
        assert region.midpoint == 1500

        region2 = GenomicRegion("chr1", 0, 501)
        assert region2.midpoint == 250  # Integer division

    def test_default_name(self):
        """Test default empty name."""
        from excavator2.io.bam import GenomicRegion

        region = GenomicRegion("chr1", 100, 200)
        assert region.name == ""


class TestReadCountResult:
    """Tests for ReadCountResult dataclass."""

    def test_creation(self):
        """Test basic ReadCountResult creation."""
        from excavator2.io.bam import GenomicRegion, ReadCountResult

        regions = [
            GenomicRegion("chr1", 1000, 2000),
            GenomicRegion("chr1", 2000, 3000),
        ]
        counts = np.array([100, 200], dtype=np.int32)

        result = ReadCountResult(
            counts=counts, regions=regions, total_reads=1000, filtered_reads=50
        )

        assert result.n_regions == 2
        assert result.total_reads == 1000
        assert result.filtered_reads == 50

    def test_mean_count(self):
        """Test mean count calculation."""
        from excavator2.io.bam import GenomicRegion, ReadCountResult

        regions = [GenomicRegion("chr1", i * 1000, (i + 1) * 1000) for i in range(3)]
        counts = np.array([100, 200, 300], dtype=np.int32)

        result = ReadCountResult(counts=counts, regions=regions, total_reads=600, filtered_reads=0)

        assert result.mean_count == 200.0

    def test_empty_result(self):
        """Test empty result."""
        from excavator2.io.bam import ReadCountResult

        result = ReadCountResult(
            counts=np.array([], dtype=np.int32), regions=[], total_reads=0, filtered_reads=0
        )

        assert result.n_regions == 0
        assert result.mean_count == 0.0  # np.mean of empty array is nan, but we handle it


class TestLoadRegionsFromBed:
    """Tests for BED file loading."""

    def test_load_simple_bed(self):
        """Test loading a simple BED file."""
        from excavator2.io.bam import load_regions_from_bed

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write("chr1\t1000\t2000\tgene1\n")
            f.write("chr1\t3000\t4000\tgene2\n")
            f.write("chr2\t1000\t2000\tgene3\n")
            bed_path = f.name

        try:
            regions = load_regions_from_bed(bed_path)

            assert len(regions) == 3
            assert regions[0].chrom == "chr1"
            assert regions[0].start == 1000
            assert regions[0].end == 2000
            assert regions[0].name == "gene1"
        finally:
            Path(bed_path).unlink()

    def test_load_bed_with_comments(self):
        """Test loading BED with comments and track lines."""
        from excavator2.io.bam import load_regions_from_bed

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write("# This is a comment\n")
            f.write("track name=test\n")
            f.write("chr1\t1000\t2000\n")
            f.write("\n")  # Empty line
            f.write("chr1\t3000\t4000\n")
            bed_path = f.name

        try:
            regions = load_regions_from_bed(bed_path)
            assert len(regions) == 2
        finally:
            Path(bed_path).unlink()

    def test_load_bed_minimal(self):
        """Test loading BED with only 3 columns."""
        from excavator2.io.bam import load_regions_from_bed

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write("chr1\t1000\t2000\n")
            bed_path = f.name

        try:
            regions = load_regions_from_bed(bed_path)
            assert len(regions) == 1
            assert regions[0].name == ""  # No name column
        finally:
            Path(bed_path).unlink()


# =============================================================================
# Tests for prepare/readcount.py - WindowData and SampleReadCounts
# =============================================================================


class TestWindowData:
    """Tests for WindowData dataclass."""

    def test_creation(self):
        """Test basic WindowData creation."""
        from excavator2.prepare.readcount import WindowData

        window = WindowData(
            chrom="chr1",
            start=1000,
            end=2000,
            position=1500,
            gc_content=0.45,
            mappability=0.95,
            region_class="IN",
            name="exon1",
        )

        assert window.chrom == "chr1"
        assert window.start == 1000
        assert window.end == 2000
        assert window.position == 1500
        assert window.gc_content == 0.45
        assert window.mappability == 0.95
        assert window.region_class == "IN"
        assert window.name == "exon1"

    def test_length_property(self):
        """Test length calculation."""
        from excavator2.prepare.readcount import WindowData

        window = WindowData(
            chrom="chr1",
            start=1000,
            end=2000,
            position=1500,
            gc_content=0.5,
            mappability=1.0,
            region_class="IN",
        )

        assert window.length == 1000


class TestSampleReadCounts:
    """Tests for SampleReadCounts dataclass."""

    def test_creation(self):
        """Test basic SampleReadCounts creation."""
        from excavator2.prepare.readcount import WindowData, SampleReadCounts

        windows = [
            WindowData("chr1", 1000, 2000, 1500, 0.5, 1.0, "IN"),
            WindowData("chr1", 2000, 3000, 2500, 0.6, 0.9, "OUT"),
            WindowData("chr2", 1000, 2000, 1500, 0.4, 1.0, "IN"),
        ]
        counts = np.array([100, 200, 150], dtype=np.int32)

        sample = SampleReadCounts(sample_name="test_sample", raw_counts=counts, windows=windows)

        assert sample.sample_name == "test_sample"
        assert sample.n_windows == 3
        assert sample.total_reads == 450
        assert set(sample.chromosomes) == {"chr1", "chr2"}

    def test_chromosome_mask(self):
        """Test chromosome mask generation."""
        from excavator2.prepare.readcount import WindowData, SampleReadCounts

        windows = [
            WindowData("chr1", 1000, 2000, 1500, 0.5, 1.0, "IN"),
            WindowData("chr1", 2000, 3000, 2500, 0.6, 0.9, "IN"),
            WindowData("chr2", 1000, 2000, 1500, 0.4, 1.0, "IN"),
        ]
        counts = np.array([100, 200, 150], dtype=np.int32)

        sample = SampleReadCounts("test", counts, windows)

        chr1_mask = sample.get_chromosome_mask("chr1")
        assert list(chr1_mask) == [True, True, False]

        chr2_mask = sample.get_chromosome_mask("chr2")
        assert list(chr2_mask) == [False, False, True]

    def test_in_out_target_masks(self):
        """Test IN/OUT target mask generation."""
        from excavator2.prepare.readcount import WindowData, SampleReadCounts

        windows = [
            WindowData("chr1", 1000, 2000, 1500, 0.5, 1.0, "IN"),
            WindowData("chr1", 2000, 3000, 2500, 0.6, 0.9, "OUT"),
            WindowData("chr1", 3000, 4000, 3500, 0.4, 1.0, "IN"),
        ]
        counts = np.array([100, 200, 150], dtype=np.int32)

        sample = SampleReadCounts("test", counts, windows)

        in_mask = sample.get_in_target_mask()
        assert list(in_mask) == [True, False, True]

        out_mask = sample.get_off_target_mask()
        assert list(out_mask) == [False, True, False]

    def test_get_chromosome_data(self):
        """Test getting data for a specific chromosome."""
        from excavator2.prepare.readcount import WindowData, SampleReadCounts

        windows = [
            WindowData("chr1", 1000, 2000, 1500, 0.5, 1.0, "IN"),
            WindowData("chr2", 1000, 2000, 1500, 0.4, 1.0, "IN"),
            WindowData("chr1", 2000, 3000, 2500, 0.6, 0.9, "IN"),
        ]
        counts = np.array([100, 200, 150], dtype=np.int32)

        sample = SampleReadCounts("test", counts, windows)

        chr1_counts, chr1_windows = sample.get_chromosome_data("chr1")
        assert len(chr1_counts) == 2
        assert list(chr1_counts) == [100, 150]
        assert len(chr1_windows) == 2


# =============================================================================
# Tests for prepare/normalize.py - Normalization algorithms
# =============================================================================


class TestMedianNormalize:
    """Tests for the median normalization function."""

    def test_basic_normalization(self):
        """Test basic median normalization."""
        from excavator2.prepare.normalize import _median_normalize

        # Create counts with bias: higher counts at higher GC
        counts = np.array([100, 110, 200, 210, 300, 310], dtype=np.float64)
        gc = np.array([20, 25, 40, 45, 60, 65], dtype=np.float64)  # GC content

        normalized, stats = _median_normalize(counts, gc, bin_size=20, min_bin_count=1)

        # After normalization, all bins should have similar median
        assert stats["n_bins"] > 0
        assert stats["master_median"] > 0
        # The variance should be reduced
        assert np.std(normalized) < np.std(counts)

    def test_empty_input(self):
        """Test normalization with empty input."""
        from excavator2.prepare.normalize import _median_normalize

        counts = np.array([], dtype=np.float64)
        feature = np.array([], dtype=np.float64)

        normalized, stats = _median_normalize(counts, feature, bin_size=5)

        assert len(normalized) == 0
        assert stats.get("n_bins", 0) == 0

    def test_all_zeros(self):
        """Test normalization when all counts are zero."""
        from excavator2.prepare.normalize import _median_normalize

        counts = np.zeros(10, dtype=np.float64)
        feature = np.arange(10, dtype=np.float64)

        normalized, stats = _median_normalize(counts, feature, bin_size=5)

        assert np.all(normalized == 0)
        assert stats["master_median"] == 0.0

    def test_uniform_data(self):
        """Test normalization with uniform data (no bias)."""
        from excavator2.prepare.normalize import _median_normalize

        np.random.seed(42)
        counts = np.random.normal(100, 10, 100).astype(np.float64)
        counts[counts < 0] = 0
        feature = np.linspace(0, 100, 100)

        normalized, stats = _median_normalize(counts, feature, bin_size=10, min_bin_count=5)

        # Should be relatively unchanged
        assert np.isclose(
            np.median(counts[counts > 0]), np.median(normalized[normalized > 0]), rtol=0.2
        )


class TestReadCountNormalizer:
    """Tests for ReadCountNormalizer class."""

    def test_normalizer_creation(self):
        """Test normalizer creation with default parameters."""
        from excavator2.prepare.normalize import ReadCountNormalizer

        normalizer = ReadCountNormalizer()

        assert normalizer.size_bin == 5.0
        assert normalizer.mappability_bin == 5.0
        assert normalizer.gc_bin == 5.0

    def test_normalizer_custom_params(self):
        """Test normalizer with custom parameters."""
        from excavator2.prepare.normalize import ReadCountNormalizer

        normalizer = ReadCountNormalizer(
            size_bin=10.0, mappability_bin=10.0, gc_bin=10.0, min_bin_count=5
        )

        assert normalizer.size_bin == 10.0
        assert normalizer.min_bin_count == 5

    def test_normalize_sample(self):
        """Test full normalization pipeline."""
        from excavator2.prepare.readcount import WindowData, SampleReadCounts
        from excavator2.prepare.normalize import ReadCountNormalizer

        # Create synthetic data with known biases
        np.random.seed(42)
        n_windows = 100

        # Simulate GC bias: higher GC -> higher counts
        gc_values = np.random.uniform(0.3, 0.7, n_windows)
        base_counts = 100 + (gc_values - 0.5) * 200  # GC bias
        counts = np.random.poisson(base_counts).astype(np.int32)

        windows = [
            WindowData(
                chrom="chr1",
                start=i * 1000,
                end=(i + 1) * 1000,
                position=i * 1000 + 500,
                gc_content=gc_values[i],
                mappability=1.0,
                region_class="IN",
            )
            for i in range(n_windows)
        ]

        sample = SampleReadCounts("test", counts, windows)
        normalizer = ReadCountNormalizer()
        result = normalizer.normalize(sample)

        # Check result structure
        assert result.sample_name == "test"
        assert result.n_windows == n_windows
        assert len(result.normalized_counts) == n_windows

        # Check that normalization reduced GC correlation
        raw_gc_corr = np.corrcoef(counts, gc_values)[0, 1]
        norm_gc_corr = np.corrcoef(result.normalized_counts, gc_values)[0, 1]
        # Correlation should be reduced (closer to 0)
        assert abs(norm_gc_corr) < abs(raw_gc_corr)

    def test_normalize_with_out_target(self):
        """Test normalization with mixed IN/OUT target regions."""
        from excavator2.prepare.readcount import WindowData, SampleReadCounts
        from excavator2.prepare.normalize import ReadCountNormalizer

        np.random.seed(42)

        windows = [
            WindowData("chr1", 0, 1000, 500, 0.5, 1.0, "IN"),
            WindowData("chr1", 1000, 2000, 1500, 0.5, 1.0, "OUT"),
            WindowData("chr1", 2000, 3000, 2500, 0.5, 1.0, "IN"),
            WindowData("chr1", 3000, 4000, 3500, 0.5, 1.0, "OUT"),
        ]
        counts = np.array([100, 200, 150, 250], dtype=np.int32)

        sample = SampleReadCounts("test", counts, windows)
        normalizer = ReadCountNormalizer()
        result = normalizer.normalize(sample)

        # Should have stats for both IN and OUT
        assert "in_target" in result.normalization_stats
        assert "off_target" in result.normalization_stats

    def test_zero_replacement(self):
        """Test that zeros are replaced with minimum non-zero value."""
        from excavator2.prepare.readcount import WindowData, SampleReadCounts
        from excavator2.prepare.normalize import ReadCountNormalizer

        windows = [
            WindowData("chr1", i * 1000, (i + 1) * 1000, i * 1000 + 500, 0.5, 1.0, "IN")
            for i in range(5)
        ]
        counts = np.array([0, 100, 200, 0, 300], dtype=np.int32)

        sample = SampleReadCounts("test", counts, windows)
        normalizer = ReadCountNormalizer()
        result = normalizer.normalize(sample)

        # No zeros should remain
        assert np.all(result.normalized_counts > 0)


class TestNormalizationIO:
    """Tests for saving/loading normalized counts."""

    def test_save_and_load_normalized_counts(self):
        """Test saving and loading normalized counts to HDF5."""
        from excavator2.prepare.readcount import WindowData
        from excavator2.prepare.normalize import (
            NormalizationResult,
            save_normalized_counts,
            load_normalized_counts,
        )

        windows = [
            WindowData("chr1", 1000, 2000, 1500, 0.5, 0.9, "IN", "gene1"),
            WindowData("chr1", 2000, 3000, 2500, 0.6, 0.8, "OUT", "gene2"),
            WindowData("chr2", 1000, 2000, 1500, 0.4, 1.0, "IN", "gene3"),
        ]
        counts = np.array([100.5, 200.3, 150.7])

        result = NormalizationResult(
            sample_name="test_sample",
            normalized_counts=counts,
            windows=windows,
            chromosomes=["chr1", "chr2"],
            normalization_stats={"test": 1.0},
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.NRC.h5"

            # Save
            save_normalized_counts(result, output_path)
            assert output_path.exists()

            # Load
            loaded = load_normalized_counts(output_path)

            assert loaded.sample_name == "test_sample"
            assert loaded.n_windows == 3
            np.testing.assert_array_almost_equal(loaded.normalized_counts, counts)
            assert len(loaded.windows) == 3
            assert loaded.windows[0].chrom == "chr1"
            assert loaded.windows[0].gc_content == 0.5


class TestReadCountIO:
    """Tests for saving/loading raw read counts."""

    def test_save_and_load_read_counts(self):
        """Test saving and loading raw read counts to HDF5."""
        from excavator2.prepare.readcount import (
            WindowData,
            SampleReadCounts,
            save_read_counts,
            load_read_counts,
        )

        windows = [
            WindowData("chr1", 1000, 2000, 1500, 0.5, 0.9, "IN", "gene1"),
            WindowData("chr2", 1000, 2000, 1500, 0.4, 1.0, "IN", "gene2"),
        ]
        counts = np.array([100, 200], dtype=np.int32)

        sample = SampleReadCounts(sample_name="test_sample", raw_counts=counts, windows=windows)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.RC.h5"

            # Save
            save_read_counts(sample, output_path)
            assert output_path.exists()

            # Load
            loaded = load_read_counts(output_path)

            assert loaded.sample_name == "test_sample"
            assert loaded.n_windows == 2
            np.testing.assert_array_equal(loaded.raw_counts, counts)


# =============================================================================
# Tests for convenience functions
# =============================================================================


class TestConvenienceFunctions:
    """Tests for module-level convenience functions."""

    def test_normalize_read_counts_function(self):
        """Test the normalize_read_counts convenience function."""
        from excavator2.prepare.readcount import WindowData, SampleReadCounts
        from excavator2.prepare import normalize_read_counts

        np.random.seed(42)
        windows = [
            WindowData("chr1", i * 1000, (i + 1) * 1000, i * 1000 + 500, 0.5, 1.0, "IN")
            for i in range(20)
        ]
        counts = np.random.poisson(100, 20).astype(np.int32)

        sample = SampleReadCounts("test", counts, windows)
        result = normalize_read_counts(sample)

        assert result.sample_name == "test"
        assert result.n_windows == 20


# =============================================================================
# Tests for module imports
# =============================================================================


class TestModuleImports:
    """Tests for module import structure."""

    def test_io_module_imports(self):
        """Test io module exports."""
        from excavator2.io import (
            BAMReader,
            GenomicRegion,
            ReadCountResult,
            count_reads_in_regions,
            load_regions_from_bed,
        )

        assert BAMReader is not None
        assert GenomicRegion is not None

    def test_prepare_module_imports(self):
        """Test prepare module exports."""
        from excavator2.prepare import (
            WindowData,
            SampleReadCounts,
            ReadCountProcessor,
            save_read_counts,
            load_read_counts,
            NormalizationResult,
            ReadCountNormalizer,
            normalize_read_counts,
            save_normalized_counts,
            load_normalized_counts,
        )

        assert WindowData is not None
        assert ReadCountNormalizer is not None
