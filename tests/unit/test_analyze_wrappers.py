"""
Tests for the Python algorithm wrappers in excavator2.analyze.

These tests verify the high-level Python API built on top of the
C++ HSLM and FastCall implementations.
"""

import pytest
import numpy as np


class TestAnalyzeImports:
    """Test that analyze module exports are available."""

    def test_import_from_analyze(self):
        """Test importing from analyze module."""
        from excavator2.analyze import (
            HSLMSegmenter,
            Segment,
            SegmentationResult,
            segment,
            FastCallCaller,
            CNVCall,
            ClassificationResult,
            CopyNumberState,
            call,
        )

    def test_import_from_excavator2(self):
        """Test importing analyze through excavator2."""
        from excavator2 import analyze
        assert hasattr(analyze, 'HSLMSegmenter')
        assert hasattr(analyze, 'FastCallCaller')


class TestHSLMSegmenter:
    """Test HSLMSegmenter wrapper class."""

    @pytest.fixture
    def segmenter(self):
        """Create default segmenter."""
        from excavator2.analyze import HSLMSegmenter
        return HSLMSegmenter()

    @pytest.fixture
    def cnv_data(self):
        """Generate data with clear CNV event."""
        np.random.seed(42)
        n = 100
        ratios = np.concatenate([
            np.random.normal(0.0, 0.1, 40),    # Normal
            np.random.normal(-0.8, 0.1, 20),   # Deletion
            np.random.normal(0.0, 0.1, 40),    # Normal
        ])
        positions = np.arange(n) * 1000
        return ratios, positions

    def test_default_parameters(self, segmenter):
        """Test default parameter values."""
        assert segmenter.omega == 0.1
        assert segmenter.theta == 1e-5
        assert segmenter.n_states == 21

    def test_custom_parameters(self):
        """Test custom parameter values."""
        from excavator2.analyze import HSLMSegmenter
        seg = HSLMSegmenter(omega=0.2, theta=1e-4, n_states=31)
        assert seg.omega == 0.2
        assert seg.theta == 1e-4
        assert seg.n_states == 31

    def test_omega_setter(self, segmenter):
        """Test omega property setter with validation."""
        segmenter.omega = 0.5
        assert segmenter.omega == 0.5

        with pytest.raises(ValueError):
            segmenter.omega = -0.1

        with pytest.raises(ValueError):
            segmenter.omega = 1.5

    def test_theta_setter(self, segmenter):
        """Test theta property setter with validation."""
        segmenter.theta = 1e-3
        assert segmenter.theta == 1e-3

        with pytest.raises(ValueError):
            segmenter.theta = 0

        with pytest.raises(ValueError):
            segmenter.theta = -1e-5

    def test_segment_returns_result(self, segmenter, cnv_data):
        """Test that segment() returns SegmentationResult."""
        from excavator2.analyze import SegmentationResult
        ratios, positions = cnv_data
        result = segmenter.segment(ratios, positions)

        assert isinstance(result, SegmentationResult)
        assert result.success
        assert result.n_segments >= 1

    def test_segment_detects_cnv(self, segmenter, cnv_data):
        """Test that segmentation detects the CNV event."""
        ratios, positions = cnv_data
        result = segmenter.segment(ratios, positions)

        assert result.success
        # Should detect at least 2 segments (normal, deletion, normal regions)
        assert result.n_segments >= 2

        # Check that segments span the full data
        total_probes = sum(s.n_probes for s in result.segments)
        assert total_probes == len(ratios)

    def test_segment_objects(self, segmenter, cnv_data):
        """Test Segment object attributes."""
        from excavator2.analyze import Segment
        ratios, positions = cnv_data
        result = segmenter.segment(ratios, positions)

        for seg in result.segments:
            assert isinstance(seg, Segment)
            assert seg.start_idx >= 0
            assert seg.end_idx > seg.start_idx
            assert seg.start_pos >= 0
            assert seg.end_pos >= seg.start_pos
            assert seg.n_probes > 0
            assert seg.n_probes == seg.end_idx - seg.start_idx

    def test_segment_with_list_input(self, segmenter):
        """Test segment() accepts Python lists."""
        np.random.seed(42)
        ratios = list(np.random.normal(0.0, 0.1, 50))
        positions = list(range(0, 50000, 1000))
        result = segmenter.segment(ratios, positions)
        assert result.success

    def test_empty_input_raises(self, segmenter):
        """Test that empty input raises ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            segmenter.segment([], [])

    def test_mismatched_lengths_raises(self, segmenter):
        """Test that mismatched input lengths raise ValueError."""
        with pytest.raises(ValueError, match="Length mismatch"):
            segmenter.segment([0.0, 0.1, 0.2], [0, 1000])

    def test_nan_input_raises(self, segmenter):
        """Test that NaN values raise ValueError."""
        with pytest.raises(ValueError, match="NaN or Inf"):
            segmenter.segment([0.0, np.nan, 0.2], [0, 1000, 2000])

    def test_inf_input_raises(self, segmenter):
        """Test that Inf values raise ValueError."""
        with pytest.raises(ValueError, match="NaN or Inf"):
            segmenter.segment([0.0, np.inf, 0.2], [0, 1000, 2000])


class TestHSLMConvenienceFunction:
    """Test segment() convenience function."""

    def test_convenience_function(self):
        """Test module-level segment() function."""
        from excavator2.analyze import segment

        np.random.seed(42)
        ratios = np.random.normal(0.0, 0.1, 50)
        positions = np.arange(50) * 1000

        result = segment(ratios, positions)
        assert result.success

    def test_convenience_function_with_params(self):
        """Test convenience function with custom parameters."""
        from excavator2.analyze import segment

        np.random.seed(42)
        ratios = np.random.normal(0.0, 0.1, 50)
        positions = np.arange(50) * 1000

        result = segment(ratios, positions, omega=0.2, theta=1e-4)
        assert result.success


class TestHSLMMultiSample:
    """Test multi-sample segmentation."""

    def test_segment_multi(self):
        """Test multi-sample segmentation."""
        from excavator2.analyze import HSLMSegmenter

        np.random.seed(42)
        # Two samples with shared breakpoint
        sample1 = np.concatenate([
            np.random.normal(0.0, 0.1, 30),
            np.random.normal(-0.5, 0.1, 20),
        ])
        sample2 = np.concatenate([
            np.random.normal(0.0, 0.1, 30),
            np.random.normal(-0.6, 0.1, 20),
        ])
        data_matrix = np.array([sample1, sample2])
        positions = np.arange(50) * 1000

        segmenter = HSLMSegmenter()
        result = segmenter.segment_multi(data_matrix, positions)

        assert result.success
        assert result.n_segments >= 1

    def test_segment_multi_invalid_dims(self):
        """Test that 1D array raises error."""
        from excavator2.analyze import HSLMSegmenter

        segmenter = HSLMSegmenter()
        with pytest.raises(ValueError, match="must be 2D"):
            segmenter.segment_multi(np.array([0.0, 0.1]), [0, 1000])


class TestFastCallCaller:
    """Test FastCallCaller wrapper class."""

    @pytest.fixture
    def caller(self):
        """Create default caller."""
        from excavator2.analyze import FastCallCaller
        return FastCallCaller()

    @pytest.fixture
    def mixed_segments(self):
        """Generate segments with mixed CN states."""
        np.random.seed(42)
        return np.concatenate([
            np.random.normal(0.0, 0.1, 20),    # Normal
            np.random.normal(-0.6, 0.1, 10),   # Deletion
            np.random.normal(0.58, 0.1, 10),   # Duplication
        ])

    def test_default_parameters(self, caller):
        """Test default parameter values."""
        assert caller.cellularity == 1.0
        assert caller.max_iterations == 1000

    def test_custom_parameters(self):
        """Test custom parameter values."""
        from excavator2.analyze import FastCallCaller
        c = FastCallCaller(cellularity=0.8, thrd=0.6, max_iterations=500)
        assert c.cellularity == 0.8
        assert c.max_iterations == 500

    def test_cellularity_validation(self):
        """Test cellularity parameter validation."""
        from excavator2.analyze import FastCallCaller

        with pytest.raises(ValueError):
            FastCallCaller(cellularity=0)

        with pytest.raises(ValueError):
            FastCallCaller(cellularity=1.5)

        with pytest.raises(ValueError):
            FastCallCaller(cellularity=-0.1)

    def test_cellularity_setter(self, caller):
        """Test cellularity property setter."""
        caller.cellularity = 0.7
        assert caller.cellularity == 0.7

        with pytest.raises(ValueError):
            caller.cellularity = 0

    def test_call_returns_result(self, caller, mixed_segments):
        """Test that call() returns ClassificationResult."""
        from excavator2.analyze import ClassificationResult
        result = caller.call(mixed_segments)

        assert isinstance(result, ClassificationResult)
        assert result.success
        assert len(result.calls) == len(mixed_segments)

    def test_call_detects_cnvs(self, caller, mixed_segments):
        """Test that classification detects CNVs."""
        result = caller.call(mixed_segments)

        assert result.success
        assert result.n_cnvs > 0

    def test_cnv_call_attributes(self, caller, mixed_segments):
        """Test CNVCall object attributes."""
        from excavator2.analyze import CNVCall, CopyNumberState
        result = caller.call(mixed_segments)

        for call in result.calls:
            assert isinstance(call, CNVCall)
            assert call.cn_call in [-2, -1, 0, 1, 2]
            assert call.absolute_cn in [0, 1, 2, 3, 4]
            assert isinstance(call.state, CopyNumberState)
            assert 0.0 <= call.probability <= 1.0

    def test_cnv_call_properties(self, caller):
        """Test CNVCall helper properties."""
        from excavator2.analyze import CNVCall, CopyNumberState

        # Normal call
        normal = CNVCall(cn_call=0, absolute_cn=2, state=CopyNumberState.NORMAL,
                        probability=0.95, segment_mean=0.0)
        assert not normal.is_cnv
        assert not normal.is_deletion
        assert not normal.is_gain

        # Deletion
        deletion = CNVCall(cn_call=-1, absolute_cn=1,
                          state=CopyNumberState.HETEROZYGOUS_DELETION,
                          probability=0.9, segment_mean=-0.5)
        assert deletion.is_cnv
        assert deletion.is_deletion
        assert not deletion.is_gain

        # Gain
        gain = CNVCall(cn_call=1, absolute_cn=3,
                      state=CopyNumberState.SINGLE_COPY_GAIN,
                      probability=0.85, segment_mean=0.5)
        assert gain.is_cnv
        assert not gain.is_deletion
        assert gain.is_gain

    def test_classification_result_helpers(self, caller, mixed_segments):
        """Test ClassificationResult helper methods."""
        result = caller.call(mixed_segments)

        # Test count properties
        assert result.n_segments == len(mixed_segments)
        assert result.n_cnvs >= 0
        assert result.n_deletions >= 0
        assert result.n_gains >= 0
        assert result.n_cnvs == result.n_deletions + result.n_gains

        # Test filter methods
        cnvs = result.get_cnvs()
        assert len(cnvs) == result.n_cnvs
        assert all(c.is_cnv for c in cnvs)

        deletions = result.get_deletions()
        assert len(deletions) == result.n_deletions
        assert all(c.is_deletion for c in deletions)

        gains = result.get_gains()
        assert len(gains) == result.n_gains
        assert all(c.is_gain for c in gains)

    def test_state_parameters(self, caller, mixed_segments):
        """Test fitted state parameters."""
        result = caller.call(mixed_segments)

        assert len(result.state_means) == 5
        assert len(result.state_sds) == 5
        assert len(result.state_priors) == 5

        # Means should be ordered
        for i in range(4):
            assert result.state_means[i] < result.state_means[i + 1]

        # SDs should be positive
        assert all(sd > 0 for sd in result.state_sds)

        # Priors should sum to ~1
        assert abs(sum(result.state_priors) - 1.0) < 0.01

    def test_em_convergence(self, caller, mixed_segments):
        """Test EM convergence information."""
        result = caller.call(mixed_segments)

        assert result.iterations > 0
        assert result.converged or result.iterations == caller.max_iterations

    def test_call_with_list_input(self, caller):
        """Test call() accepts Python lists."""
        np.random.seed(42)
        segment_means = list(np.random.normal(0.0, 0.1, 30))
        result = caller.call(segment_means)
        assert result.success

    def test_empty_input_raises(self, caller):
        """Test that empty input raises ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            caller.call([])

    def test_nan_input_raises(self, caller):
        """Test that NaN values raise ValueError."""
        with pytest.raises(ValueError, match="NaN or Inf"):
            caller.call([0.0, np.nan, 0.2])


class TestFastCallConvenienceFunction:
    """Test call() convenience function."""

    def test_convenience_function(self):
        """Test module-level call() function."""
        from excavator2.analyze import call

        np.random.seed(42)
        segment_means = np.random.normal(0.0, 0.1, 30)
        result = call(segment_means)
        assert result.success

    def test_convenience_function_with_params(self):
        """Test convenience function with custom parameters."""
        from excavator2.analyze import call

        np.random.seed(42)
        segment_means = np.random.normal(0.0, 0.1, 30)
        result = call(segment_means, cellularity=0.8, thrd=0.6)
        assert result.success


class TestCopyNumberState:
    """Test CopyNumberState enum."""

    def test_enum_values(self):
        """Test enum values are correct."""
        from excavator2.analyze import CopyNumberState

        assert CopyNumberState.HOMOZYGOUS_DELETION.value == -2
        assert CopyNumberState.HETEROZYGOUS_DELETION.value == -1
        assert CopyNumberState.NORMAL.value == 0
        assert CopyNumberState.SINGLE_COPY_GAIN.value == 1
        assert CopyNumberState.AMPLIFICATION.value == 2

    def test_from_absolute_cn(self):
        """Test converting absolute CN to state."""
        from excavator2.analyze import CopyNumberState

        assert CopyNumberState.from_absolute_cn(0) == CopyNumberState.HOMOZYGOUS_DELETION
        assert CopyNumberState.from_absolute_cn(1) == CopyNumberState.HETEROZYGOUS_DELETION
        assert CopyNumberState.from_absolute_cn(2) == CopyNumberState.NORMAL
        assert CopyNumberState.from_absolute_cn(3) == CopyNumberState.SINGLE_COPY_GAIN
        assert CopyNumberState.from_absolute_cn(4) == CopyNumberState.AMPLIFICATION

    def test_absolute_cn_property(self):
        """Test absolute_cn property."""
        from excavator2.analyze import CopyNumberState

        assert CopyNumberState.HOMOZYGOUS_DELETION.absolute_cn == 0
        assert CopyNumberState.NORMAL.absolute_cn == 2
        assert CopyNumberState.AMPLIFICATION.absolute_cn == 4

    def test_label_property(self):
        """Test label property."""
        from excavator2.analyze import CopyNumberState

        assert "CN=0" in CopyNumberState.HOMOZYGOUS_DELETION.label
        assert "CN=2" in CopyNumberState.NORMAL.label
        assert "CN=4" in CopyNumberState.AMPLIFICATION.label


class TestEndToEndPipeline:
    """Test full segmentation + calling pipeline."""

    def test_pipeline_integration(self):
        """Test HSLM segmentation followed by FastCall classification."""
        from excavator2.analyze import HSLMSegmenter, FastCallCaller

        np.random.seed(42)

        # Generate data with CNV event
        n = 100
        ratios = np.concatenate([
            np.random.normal(0.0, 0.08, 40),    # Normal
            np.random.normal(-0.8, 0.08, 20),   # Deletion
            np.random.normal(0.0, 0.08, 40),    # Normal
        ])
        positions = np.arange(n) * 1000

        # Step 1: Segment
        segmenter = HSLMSegmenter(omega=0.1, theta=1e-5)
        seg_result = segmenter.segment(ratios, positions)

        assert seg_result.success
        assert seg_result.n_segments >= 2

        # Step 2: Call CNVs on segment means
        segment_means = [s.mean for s in seg_result.segments]

        caller = FastCallCaller(cellularity=1.0)
        call_result = caller.call(segment_means)

        assert call_result.success
        assert len(call_result.calls) == len(seg_result.segments)

        # Should detect at least one deletion
        assert call_result.n_deletions >= 1

    def test_pipeline_all_normal(self):
        """Test pipeline on data with no CNVs."""
        from excavator2.analyze import HSLMSegmenter, FastCallCaller

        np.random.seed(123)
        ratios = np.random.normal(0.0, 0.1, 50)
        positions = np.arange(50) * 1000

        segmenter = HSLMSegmenter()
        seg_result = segmenter.segment(ratios, positions)
        assert seg_result.success

        segment_means = [s.mean for s in seg_result.segments]
        caller = FastCallCaller()
        call_result = caller.call(segment_means)

        assert call_result.success
        # Most or all segments should be normal
        normal_count = sum(1 for c in call_result.calls if c.cn_call == 0)
        assert normal_count >= len(call_result.calls) * 0.8
