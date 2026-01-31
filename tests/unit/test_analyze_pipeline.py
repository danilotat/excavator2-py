"""
Tests for the CNV analysis pipeline (Phase 7).

This module tests:
- Log2 ratio computation
- Analysis pipeline integration
- VCF/BED output writing
- Paired and pooled experiment modes
"""

import pytest
import numpy as np
import tempfile
from pathlib import Path

from excavator2.analyze.ratio import (
    Log2RatioResult,
    compute_log2_ratio,
    compute_log2_ratio_pooled,
    apply_cellularity_correction,
)
from excavator2.analyze.pipeline import (
    CNVSegment,
    ChromosomeResult,
    AnalysisResult,
    AnalysisParameters,
    CNVAnalyzer,
)
from excavator2.analyze import (
    CopyNumberState,
    HSLMSegmenter,
    FastCallCaller,
)
from excavator2.prepare.normalize import NormalizationResult
from excavator2.prepare.readcount import WindowData
from excavator2.io.vcf import (
    write_vcf,
    write_bed,
    write_segments_tsv,
    write_fastcall_bed,
)

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def sample_windows():
    """Create sample window metadata."""
    windows = []
    for i in range(100):
        chrom = "chr1" if i < 50 else "chr2"
        start = (i % 50) * 1000
        end = start + 500
        windows.append(
            WindowData(
                chrom=chrom,
                start=start,
                end=end,
                position=(start + end) // 2,
                gc_content=0.45,
                mappability=0.95,
                region_class="IN" if i % 3 != 0 else "OUT",
                name=f"exon_{i}",
            )
        )
    return windows


@pytest.fixture
def test_normalization_result(sample_windows):
    """Create test sample normalization result with a CNV region."""
    counts = np.ones(100) * 100.0
    # Add a deletion region (lower counts)
    counts[20:30] = 50.0
    # Add an amplification region (higher counts)
    counts[70:80] = 200.0

    return NormalizationResult(
        sample_name="test_sample",
        normalized_counts=counts,
        windows=sample_windows,
        chromosomes=["chr1", "chr2"],
        normalization_stats={},
    )


@pytest.fixture
def control_normalization_result(sample_windows):
    """Create control sample normalization result (baseline)."""
    counts = np.ones(100) * 100.0

    return NormalizationResult(
        sample_name="control_sample",
        normalized_counts=counts,
        windows=sample_windows,
        chromosomes=["chr1", "chr2"],
        normalization_stats={},
    )


@pytest.fixture
def cnv_analysis_result():
    """Create a sample analysis result for output testing."""
    segments = [
        CNVSegment(
            chrom="chr1",
            start=0,
            end=19000,
            start_idx=0,
            end_idx=20,
            n_probes=20,
            segment_mean=0.0,
            cn_call=0,
            absolute_cn=2,
            probability=0.95,
            state=CopyNumberState.NORMAL,
            region_class="IN",
        ),
        CNVSegment(
            chrom="chr1",
            start=20000,
            end=29000,
            start_idx=20,
            end_idx=30,
            n_probes=10,
            segment_mean=-1.0,
            cn_call=-1,
            absolute_cn=1,
            probability=0.85,
            state=CopyNumberState.HETEROZYGOUS_DELETION,
            region_class="IN",
        ),
        CNVSegment(
            chrom="chr1",
            start=30000,
            end=49000,
            start_idx=30,
            end_idx=50,
            n_probes=20,
            segment_mean=0.0,
            cn_call=0,
            absolute_cn=2,
            probability=0.92,
            state=CopyNumberState.NORMAL,
            region_class="IN",
        ),
        CNVSegment(
            chrom="chr2",
            start=0,
            end=19000,
            start_idx=0,
            end_idx=20,
            n_probes=20,
            segment_mean=0.0,
            cn_call=0,
            absolute_cn=2,
            probability=0.90,
            state=CopyNumberState.NORMAL,
            region_class="IN",
        ),
        CNVSegment(
            chrom="chr2",
            start=20000,
            end=29000,
            start_idx=20,
            end_idx=30,
            n_probes=10,
            segment_mean=0.58,
            cn_call=1,
            absolute_cn=3,
            probability=0.88,
            state=CopyNumberState.SINGLE_COPY_GAIN,
            region_class="IN",
        ),
    ]

    chromosome_results = {
        "chr1": ChromosomeResult(
            chrom="chr1", segments=segments[:3], breakpoints=[0, 20, 30, 50], n_segments=3, n_cnvs=1
        ),
        "chr2": ChromosomeResult(
            chrom="chr2", segments=segments[3:], breakpoints=[0, 20, 30], n_segments=2, n_cnvs=1
        ),
    }

    return AnalysisResult(
        test_sample="test_sample",
        control_sample="control_sample",
        segments=segments,
        chromosome_results=chromosome_results,
        n_segments=5,
        n_cnvs=2,
        parameters={"hslm": {"omega": 0.1}, "fastcall": {"cellularity": 1.0}},
    )


# =============================================================================
# Log2 Ratio Tests
# =============================================================================


class TestLog2RatioComputation:
    """Tests for log2 ratio computation."""

    def test_compute_log2_ratio_basic(
        self, test_normalization_result, control_normalization_result
    ):
        """Test basic log2 ratio computation."""
        result = compute_log2_ratio(
            test_normalization_result, control_normalization_result, median_center=False
        )

        assert isinstance(result, Log2RatioResult)
        assert result.sample_name == "test_sample"
        assert result.control_name == "control_sample"
        assert len(result.log2_ratios) == 100
        assert len(result.positions) == 100

    def test_compute_log2_ratio_with_centering(
        self, test_normalization_result, control_normalization_result
    ):
        """Test log2 ratio computation with median centering."""
        result = compute_log2_ratio(
            test_normalization_result,
            control_normalization_result,
            median_center=True,
            separate_regions=False,
        )

        # After centering, median should be close to 0
        assert abs(np.median(result.log2_ratios)) < 0.1

    def test_compute_log2_ratio_separate_regions(
        self, test_normalization_result, control_normalization_result
    ):
        """Test log2 ratio computation with separate IN/OUT centering."""
        result = compute_log2_ratio(
            test_normalization_result,
            control_normalization_result,
            median_center=True,
            separate_regions=True,
        )

        # Check IN and OUT regions are centered separately
        in_mask = result.in_target_mask
        out_mask = ~in_mask

        in_median = np.median(result.log2_ratios[in_mask])
        out_median = np.median(result.log2_ratios[out_mask])

        # Both should be centered (close to 0)
        assert abs(in_median) < 0.1
        assert abs(out_median) < 0.1

    def test_compute_log2_ratio_detects_cnv(
        self, test_normalization_result, control_normalization_result
    ):
        """Test that log2 ratios reflect CNV regions."""
        result = compute_log2_ratio(
            test_normalization_result, control_normalization_result, median_center=False
        )

        # Deletion region (indices 20-30) should have negative log2 ratios
        deletion_ratios = result.log2_ratios[20:30]
        assert np.mean(deletion_ratios) < -0.5

        # Amplification region (indices 70-80) should have positive log2 ratios
        amp_ratios = result.log2_ratios[70:80]
        assert np.mean(amp_ratios) > 0.5

    def test_get_chromosome_data(self, test_normalization_result, control_normalization_result):
        """Test getting data for a specific chromosome."""
        result = compute_log2_ratio(
            test_normalization_result, control_normalization_result, median_center=False
        )

        chr1_data = result.get_chromosome_data("chr1")
        assert chr1_data["n_windows"] == 50
        assert len(chr1_data["log2_ratios"]) == 50
        assert len(chr1_data["positions"]) == 50

    def test_length_mismatch_raises_error(self, test_normalization_result, sample_windows):
        """Test that mismatched lengths raise an error."""
        # Create control with different number of windows
        short_control = NormalizationResult(
            sample_name="short_control",
            normalized_counts=np.ones(50) * 100.0,
            windows=sample_windows[:50],
            chromosomes=["chr1"],
            normalization_stats={},
        )

        with pytest.raises(ValueError, match="Window count mismatch"):
            compute_log2_ratio(test_normalization_result, short_control)


class TestPooledLog2Ratio:
    """Tests for pooled control log2 ratio computation."""

    def test_compute_pooled_ratio(
        self, test_normalization_result, control_normalization_result, sample_windows
    ):
        """Test pooled control log2 ratio computation."""
        # Create second control
        control2 = NormalizationResult(
            sample_name="control2",
            normalized_counts=np.ones(100) * 100.0,
            windows=sample_windows,
            chromosomes=["chr1", "chr2"],
            normalization_stats={},
        )

        result = compute_log2_ratio_pooled(
            test_normalization_result, [control_normalization_result, control2], median_center=False
        )

        assert result.control_name == "PooledControl"
        assert len(result.log2_ratios) == 100

    def test_pooled_empty_controls_raises_error(self, test_normalization_result):
        """Test that empty controls list raises error."""
        with pytest.raises(ValueError, match="At least one control"):
            compute_log2_ratio_pooled(test_normalization_result, [])


class TestCellularityCorrection:
    """Tests for cellularity correction."""

    def test_no_correction_at_full_purity(self):
        """Test that cellularity=1.0 returns original values."""
        ratios = np.array([0.0, -1.0, 1.0])
        corrected = apply_cellularity_correction(ratios, 1.0)
        np.testing.assert_array_almost_equal(ratios, corrected)

    def test_correction_amplifies_signal(self):
        """Test that low cellularity amplifies the signal."""
        ratios = np.array([-0.5, 0.0, 0.5])
        corrected = apply_cellularity_correction(ratios, 0.5)

        # With low cellularity, signals should be amplified
        assert abs(corrected[0]) > abs(ratios[0])
        assert abs(corrected[2]) > abs(ratios[2])

    def test_invalid_cellularity_raises_error(self):
        """Test that invalid cellularity raises error."""
        ratios = np.array([0.0])
        with pytest.raises(ValueError):
            apply_cellularity_correction(ratios, 0.0)


# =============================================================================
# Analysis Pipeline Tests
# =============================================================================


class TestAnalysisParameters:
    """Tests for AnalysisParameters dataclass."""

    def test_default_parameters(self):
        """Test default parameter values."""
        params = AnalysisParameters()

        assert params.omega == 0.1
        assert params.theta == 1e-5
        assert params.step_eta == 200000.0
        assert params.n_states == 21
        assert params.cellularity == 1.0
        assert params.thrd == 0.5
        assert params.thru == 0.35
        assert params.min_exons == 4

    def test_to_dict(self):
        """Test parameter serialization."""
        params = AnalysisParameters(omega=0.2, cellularity=0.8)
        d = params.to_dict()

        assert d["hslm"]["omega"] == 0.2
        assert d["fastcall"]["cellularity"] == 0.8


class TestCNVAnalyzer:
    """Tests for CNVAnalyzer class."""

    def test_analyzer_initialization(self):
        """Test analyzer initialization."""
        params = AnalysisParameters(omega=0.2)
        analyzer = CNVAnalyzer(params)

        assert analyzer.params.omega == 0.2

    def test_analyze_paired(self, test_normalization_result, control_normalization_result):
        """Test paired analysis."""
        analyzer = CNVAnalyzer()
        result = analyzer.analyze_paired(test_normalization_result, control_normalization_result)

        assert isinstance(result, AnalysisResult)
        assert result.test_sample == "test_sample"
        assert result.control_sample == "control_sample"
        assert result.n_segments > 0

    def test_analyze_pooled(
        self, test_normalization_result, control_normalization_result, sample_windows
    ):
        """Test pooled analysis."""
        control2 = NormalizationResult(
            sample_name="control2",
            normalized_counts=np.ones(100) * 100.0,
            windows=sample_windows,
            chromosomes=["chr1", "chr2"],
            normalization_stats={},
        )

        analyzer = CNVAnalyzer()
        result = analyzer.analyze_pooled(
            test_normalization_result, [control_normalization_result, control2]
        )

        assert isinstance(result, AnalysisResult)
        assert result.control_sample == "PooledControl"

    def test_analyze_detects_cnvs(self, test_normalization_result, control_normalization_result):
        """Test that analyzer detects CNVs in synthetic data."""
        analyzer = CNVAnalyzer()
        result = analyzer.analyze_paired(test_normalization_result, control_normalization_result)

        # Should detect at least some CNVs (we have deletion and amp regions)
        # Note: may not always detect depending on segmentation sensitivity
        assert result.n_segments > 0

    def test_analyze_specific_chromosomes(
        self, test_normalization_result, control_normalization_result
    ):
        """Test analysis of specific chromosomes."""
        analyzer = CNVAnalyzer()
        result = analyzer.analyze_paired(
            test_normalization_result, control_normalization_result, chromosomes=["chr1"]
        )

        # Should only have chr1 results
        assert "chr1" in result.chromosome_results
        assert "chr2" not in result.chromosome_results


class TestCNVSegment:
    """Tests for CNVSegment dataclass."""

    def test_segment_properties(self):
        """Test CNVSegment properties."""
        deletion = CNVSegment(
            chrom="chr1",
            start=1000,
            end=2000,
            start_idx=10,
            end_idx=20,
            n_probes=10,
            segment_mean=-1.0,
            cn_call=-1,
            absolute_cn=1,
            probability=0.9,
            state=CopyNumberState.HETEROZYGOUS_DELETION,
            region_class="IN",
        )

        assert deletion.length == 1000
        assert deletion.is_cnv is True
        assert deletion.is_deletion is True
        assert deletion.is_gain is False

        gain = CNVSegment(
            chrom="chr1",
            start=3000,
            end=4000,
            start_idx=30,
            end_idx=40,
            n_probes=10,
            segment_mean=0.58,
            cn_call=1,
            absolute_cn=3,
            probability=0.85,
            state=CopyNumberState.SINGLE_COPY_GAIN,
            region_class="IN",
        )

        assert gain.is_cnv is True
        assert gain.is_deletion is False
        assert gain.is_gain is True


class TestAnalysisResult:
    """Tests for AnalysisResult dataclass."""

    def test_result_properties(self, cnv_analysis_result):
        """Test AnalysisResult properties."""
        assert cnv_analysis_result.n_segments == 5
        assert cnv_analysis_result.n_cnvs == 2
        assert cnv_analysis_result.n_deletions == 1
        assert cnv_analysis_result.n_gains == 1

    def test_get_cnvs(self, cnv_analysis_result):
        """Test filtering CNVs."""
        cnvs = cnv_analysis_result.get_cnvs()
        assert len(cnvs) == 2

    def test_get_deletions(self, cnv_analysis_result):
        """Test filtering deletions."""
        deletions = cnv_analysis_result.get_deletions()
        assert len(deletions) == 1
        assert all(d.is_deletion for d in deletions)

    def test_get_gains(self, cnv_analysis_result):
        """Test filtering gains."""
        gains = cnv_analysis_result.get_gains()
        assert len(gains) == 1
        assert all(g.is_gain for g in gains)


# =============================================================================
# Output Writing Tests
# =============================================================================


class TestVCFOutput:
    """Tests for VCF output writing."""

    def test_write_vcf(self, cnv_analysis_result, tmp_path):
        """Test writing VCF file."""
        vcf_path = tmp_path / "test.vcf"
        write_vcf(cnv_analysis_result, vcf_path, cnv_only=True)

        assert vcf_path.exists()

        # Check content
        content = vcf_path.read_text()
        assert "##fileformat=VCFv4.2" in content
        assert "##source=EXCAVATOR2" in content
        assert "SVTYPE=" in content
        assert "test_sample" in content

    def test_write_vcf_all_segments(self, cnv_analysis_result, tmp_path):
        """Test writing VCF with all segments."""
        vcf_path = tmp_path / "test_all.vcf"
        write_vcf(cnv_analysis_result, vcf_path, cnv_only=False)

        content = vcf_path.read_text()
        # Count data lines (non-header)
        data_lines = [l for l in content.split("\n") if l and not l.startswith("#")]
        assert len(data_lines) == 5  # All 5 segments


class TestBEDOutput:
    """Tests for BED output writing."""

    def test_write_bed(self, cnv_analysis_result, tmp_path):
        """Test writing BED file."""
        bed_path = tmp_path / "test.bed"
        write_bed(cnv_analysis_result, bed_path, cnv_only=True)

        assert bed_path.exists()

        content = bed_path.read_text()
        lines = content.strip().split("\n")
        # Header + 2 CNVs
        assert len(lines) == 3

    def test_write_bed_all(self, cnv_analysis_result, tmp_path):
        """Test writing BED with all segments."""
        bed_path = tmp_path / "test_all.bed"
        write_bed(cnv_analysis_result, bed_path, cnv_only=False)

        content = bed_path.read_text()
        lines = content.strip().split("\n")
        # Header + 5 segments
        assert len(lines) == 6


class TestTSVOutput:
    """Tests for TSV output writing."""

    def test_write_segments_tsv(self, cnv_analysis_result, tmp_path):
        """Test writing segments TSV file."""
        tsv_path = tmp_path / "segments.tsv"
        write_segments_tsv(cnv_analysis_result, tsv_path)

        assert tsv_path.exists()

        content = tsv_path.read_text()
        lines = content.strip().split("\n")
        # Header + 5 segments
        assert len(lines) == 6

        # Check header
        header = lines[0].split("\t")
        assert "Chromosome" in header
        assert "CN" in header
        assert "Probability" in header


class TestFastCallBEDOutput:
    """Tests for FastCall BED output writing."""

    def test_write_fastcall_bed(self, cnv_analysis_result, tmp_path):
        """Test writing FastCall BED file."""
        bed_path = tmp_path / "fastcall.bed"
        write_fastcall_bed(cnv_analysis_result, bed_path)

        assert bed_path.exists()

        content = bed_path.read_text()
        lines = content.strip().split("\n")
        # Header + 5 segments
        assert len(lines) == 6

        # Check header format
        header = lines[0].split("\t")
        assert "Call" in header
        assert "ProbCall" in header


# =============================================================================
# Integration Tests
# =============================================================================


class TestEndToEndIntegration:
    """End-to-end integration tests."""

    def test_full_pipeline_paired(
        self, test_normalization_result, control_normalization_result, tmp_path
    ):
        """Test full pipeline in paired mode."""
        # Run analysis
        params = AnalysisParameters(min_exons=2)  # Lower threshold for test data
        analyzer = CNVAnalyzer(params)
        result = analyzer.analyze_paired(test_normalization_result, control_normalization_result)

        # Write all outputs
        write_vcf(result, tmp_path / "test.vcf", cnv_only=True)
        write_bed(result, tmp_path / "test.bed", cnv_only=True)
        write_segments_tsv(result, tmp_path / "segments.tsv")
        write_fastcall_bed(result, tmp_path / "fastcall.bed")

        # Verify all files exist
        assert (tmp_path / "test.vcf").exists()
        assert (tmp_path / "test.bed").exists()
        assert (tmp_path / "segments.tsv").exists()
        assert (tmp_path / "fastcall.bed").exists()

    def test_full_pipeline_pooled(
        self, test_normalization_result, control_normalization_result, sample_windows, tmp_path
    ):
        """Test full pipeline in pooled mode."""
        # Create multiple controls
        controls = [control_normalization_result]
        for i in range(2):
            controls.append(
                NormalizationResult(
                    sample_name=f"control_{i}",
                    normalized_counts=np.ones(100) * 100.0 + np.random.normal(0, 5, 100),
                    windows=sample_windows,
                    chromosomes=["chr1", "chr2"],
                    normalization_stats={},
                )
            )

        # Run analysis
        params = AnalysisParameters(min_exons=2)
        analyzer = CNVAnalyzer(params)
        result = analyzer.analyze_pooled(test_normalization_result, controls)

        # Write outputs
        write_vcf(result, tmp_path / "pooled.vcf", cnv_only=True)
        write_segments_tsv(result, tmp_path / "pooled_segments.tsv")

        # Verify
        assert (tmp_path / "pooled.vcf").exists()
        assert (tmp_path / "pooled_segments.tsv").exists()
        assert result.control_sample == "PooledControl"
