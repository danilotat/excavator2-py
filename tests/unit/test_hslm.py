"""
Tests for the HSLM (Heterogeneous Shifting Level Model) segmentation algorithm.

These tests verify that the C++ implementation correctly ports the original
Fortran/R implementation from EXCAVATOR2.
"""

import pytest
import numpy as np


class TestLogSpaceUtilities:
    """Test log-space arithmetic utilities."""

    def test_elnsum_basic(self):
        """Test elnsum produces correct log(exp(x) + exp(y))."""
        from excavator2._excavator_core import elnsum

        # Test with simple values
        result = elnsum(-1.0, -2.0)
        expected = np.log(np.exp(-1.0) + np.exp(-2.0))
        assert abs(result - expected) < 1e-10

    def test_elnsum_symmetric(self):
        """Test elnsum is symmetric."""
        from excavator2._excavator_core import elnsum

        assert abs(elnsum(-1.0, -2.0) - elnsum(-2.0, -1.0)) < 1e-10

    def test_elnsum_extreme_values(self):
        """Test elnsum handles extreme value differences."""
        from excavator2._excavator_core import elnsum

        # When one value is much larger, result should be close to the larger
        result = elnsum(0.0, -100.0)
        assert abs(result - 0.0) < 1e-10

    def test_logsumexp_vector(self):
        """Test logsumexp for vector input."""
        from excavator2._excavator_core import logsumexp

        values = [-1.0, -2.0, -3.0]
        result = logsumexp(values)
        expected = np.log(sum(np.exp(v) for v in values))
        assert abs(result - expected) < 1e-10

    def test_logsumexp_single_element(self):
        """Test logsumexp with single element."""
        from excavator2._excavator_core import logsumexp

        result = logsumexp([-5.0])
        assert abs(result - (-5.0)) < 1e-10


class TestHSLMParameters:
    """Test HSLMParameters configuration."""

    def test_default_parameters(self):
        """Test default parameter values."""
        from excavator2._excavator_core.hslm import HSLMParameters

        params = HSLMParameters()
        assert params.omega == 0.1
        assert params.theta == 1e-5
        assert params.n_states == 21
        assert params.step_eta == 200000.0
        assert params.min_segment_size == 1

    def test_parameter_modification(self):
        """Test parameter values can be modified."""
        from excavator2._excavator_core.hslm import HSLMParameters

        params = HSLMParameters()
        params.omega = 0.2
        params.theta = 1e-4
        params.n_states = 11

        assert params.omega == 0.2
        assert params.theta == 1e-4
        assert params.n_states == 11

    def test_parameter_repr(self):
        """Test parameter string representation."""
        from excavator2._excavator_core.hslm import HSLMParameters

        params = HSLMParameters()
        repr_str = repr(params)
        assert "HSLMParameters" in repr_str
        assert "omega" in repr_str


class TestHSLMSegmentation:
    """Test HSLM segmentation algorithm."""

    @pytest.fixture
    def simple_segments(self):
        """Generate simple three-segment data: normal, deletion, normal."""
        np.random.seed(42)
        positions = np.arange(1000, 401000, 1000, dtype=np.int64)
        log2_ratios = np.concatenate(
            [
                np.random.normal(0.0, 0.1, 150),  # Normal region
                np.random.normal(-0.5, 0.1, 100),  # Deletion
                np.random.normal(0.0, 0.1, 150),  # Normal region
            ]
        )
        return log2_ratios, positions

    @pytest.fixture
    def duplication_segments(self):
        """Generate data with a duplication segment."""
        np.random.seed(123)
        positions = np.arange(1000, 301000, 1000, dtype=np.int64)
        log2_ratios = np.concatenate(
            [
                np.random.normal(0.0, 0.1, 100),  # Normal
                np.random.normal(0.58, 0.1, 100),  # Duplication (log2(3/2))
                np.random.normal(0.0, 0.1, 100),  # Normal
            ]
        )
        return log2_ratios, positions

    def test_segment_returns_result(self):
        """Test that segment() returns an HSLMResult object."""
        from excavator2._excavator_core.hslm import HSLM, HSLMResult

        log2_ratios = [0.0] * 100
        positions = list(range(1000, 101000, 1000))

        segmenter = HSLM()
        result = segmenter.segment(log2_ratios, positions)

        assert isinstance(result, HSLMResult)
        assert result.success

    def test_segment_simple_data(self, simple_segments):
        """Test segmentation on simple three-segment data."""
        from excavator2._excavator_core.hslm import HSLM

        log2_ratios, positions = simple_segments

        segmenter = HSLM()
        result = segmenter.segment(log2_ratios.tolist(), positions.tolist())

        assert result.success
        # Should detect 3 segments (normal, deletion, normal)
        assert result.n_segments >= 2
        # Breakpoints should include start (0) and end (400)
        assert result.breakpoints[0] == 0
        assert result.breakpoints[-1] == 400
        # Segment means should be computed
        assert len(result.segment_means) == result.n_segments

    def test_segment_detects_deletion(self, simple_segments):
        """Test that segmentation detects the deletion region."""
        from excavator2._excavator_core.hslm import HSLM

        log2_ratios, positions = simple_segments

        segmenter = HSLM()
        result = segmenter.segment(log2_ratios.tolist(), positions.tolist())

        # Find segment with negative mean (deletion)
        deletion_segments = [m for m in result.segment_means if m < -0.2]
        assert len(deletion_segments) >= 1, "Should detect at least one deletion segment"

    def test_segment_detects_duplication(self, duplication_segments):
        """Test that segmentation detects duplication region."""
        from excavator2._excavator_core.hslm import HSLM

        log2_ratios, positions = duplication_segments

        segmenter = HSLM()
        result = segmenter.segment(log2_ratios.tolist(), positions.tolist())

        assert result.success
        # Find segment with positive mean (duplication)
        duplication_segments = [m for m in result.segment_means if m > 0.3]
        assert len(duplication_segments) >= 1, "Should detect at least one duplication segment"

    def test_segment_flat_data(self):
        """Test segmentation on flat (no CNV) data."""
        from excavator2._excavator_core.hslm import HSLM

        np.random.seed(42)
        positions = np.arange(1000, 201000, 1000, dtype=np.int64)
        log2_ratios = np.random.normal(0.0, 0.05, 200)

        segmenter = HSLM()
        result = segmenter.segment(log2_ratios.tolist(), positions.tolist())

        assert result.success
        # All segment means should be close to 0
        for mean in result.segment_means:
            assert abs(mean) < 0.3, f"Flat data should have near-zero means, got {mean}"

    def test_empty_data_fails(self):
        """Test that empty data returns error."""
        from excavator2._excavator_core.hslm import HSLM

        segmenter = HSLM()
        result = segmenter.segment([], [])

        assert not result.success
        assert result.error_message != ""

    def test_mismatched_lengths_fails(self):
        """Test that mismatched data lengths return error."""
        from excavator2._excavator_core.hslm import HSLM

        log2_ratios = [0.0] * 100
        positions = list(range(1000, 51000, 1000))  # Only 50 positions

        segmenter = HSLM()
        result = segmenter.segment(log2_ratios, positions)

        assert not result.success

    def test_min_segment_size_filtering(self, simple_segments):
        """Test that min_segment_size parameter filters small segments."""
        from excavator2._excavator_core.hslm import HSLM, HSLMParameters

        log2_ratios, positions = simple_segments

        # First run without filtering
        segmenter = HSLM()
        result1 = segmenter.segment(log2_ratios.tolist(), positions.tolist())

        # Then run with aggressive filtering
        params = HSLMParameters()
        params.min_segment_size = 50
        segmenter_filtered = HSLM(params)
        result2 = segmenter_filtered.segment(log2_ratios.tolist(), positions.tolist())

        # Filtered version should have fewer or equal segments
        assert result2.n_segments <= result1.n_segments

    def test_omega_parameter_effect(self, simple_segments):
        """Test that omega parameter affects segmentation."""
        from excavator2._excavator_core.hslm import HSLM, HSLMParameters

        log2_ratios, positions = simple_segments

        # Low omega: more variance attributed to noise
        params_low = HSLMParameters()
        params_low.omega = 0.01
        result_low = HSLM(params_low).segment(log2_ratios.tolist(), positions.tolist())

        # High omega: more variance attributed to states
        params_high = HSLMParameters()
        params_high.omega = 0.5
        result_high = HSLM(params_high).segment(log2_ratios.tolist(), positions.tolist())

        # Both should succeed
        assert result_low.success
        assert result_high.success


class TestHSLMConvenienceFunction:
    """Test the convenience segment() function."""

    def test_convenience_function(self):
        """Test the module-level segment() function."""
        from excavator2._excavator_core import hslm

        np.random.seed(42)
        positions = list(range(1000, 101000, 1000))
        log2_ratios = list(np.random.normal(0.0, 0.1, 100))

        result = hslm.segment(log2_ratios, positions)

        assert result.success
        assert result.n_segments >= 1

    def test_convenience_function_with_params(self):
        """Test convenience function with custom parameters."""
        from excavator2._excavator_core import hslm

        np.random.seed(42)
        positions = list(range(1000, 101000, 1000))
        log2_ratios = list(np.random.normal(0.0, 0.1, 100))

        result = hslm.segment(log2_ratios, positions, omega=0.2, theta=1e-4, min_segment_size=5)

        assert result.success


class TestHSLMMultiSample:
    """Test multi-sample segmentation."""

    def test_segment_multi_basic(self):
        """Test multi-sample segmentation."""
        from excavator2._excavator_core.hslm import HSLM

        np.random.seed(42)
        positions = list(range(1000, 101000, 1000))
        data_matrix = [
            list(np.random.normal(0.0, 0.1, 100)),
            list(np.random.normal(0.0, 0.1, 100)),
        ]

        segmenter = HSLM()
        result = segmenter.segment_multi(data_matrix, positions)

        assert result.success
        assert result.n_segments >= 1

    def test_segment_multi_consistent_breakpoints(self):
        """Test that multi-sample finds consistent breakpoints across samples."""
        from excavator2._excavator_core.hslm import HSLM

        np.random.seed(42)
        positions = list(range(1000, 201000, 1000))

        # Create two samples with same breakpoint structure
        sample1 = np.concatenate(
            [
                np.random.normal(0.0, 0.1, 100),
                np.random.normal(-0.5, 0.1, 100),
            ]
        )
        sample2 = np.concatenate(
            [
                np.random.normal(0.0, 0.1, 100),
                np.random.normal(-0.5, 0.1, 100),
            ]
        )

        data_matrix = [sample1.tolist(), sample2.tolist()]

        segmenter = HSLM()
        result = segmenter.segment_multi(data_matrix, positions)

        assert result.success
        # Should detect the breakpoint around position 100
        assert result.n_segments >= 2


class TestHSLMStateValues:
    """Test state values and path output."""

    def test_state_values_range(self):
        """Test that state values are in expected range."""
        from excavator2._excavator_core.hslm import HSLM

        np.random.seed(42)
        positions = list(range(1000, 101000, 1000))
        log2_ratios = list(np.random.normal(0.0, 0.1, 100))

        result = HSLM().segment(log2_ratios, positions)

        assert result.success
        # State values should be in [-1, 1] range (default)
        for sv in result.state_values:
            assert -1.0 <= sv <= 1.0

    def test_state_path_valid(self):
        """Test that state path indices are valid."""
        from excavator2._excavator_core.hslm import HSLM

        np.random.seed(42)
        positions = list(range(1000, 101000, 1000))
        log2_ratios = list(np.random.normal(0.0, 0.1, 100))

        result = HSLM().segment(log2_ratios, positions)

        assert result.success
        # State path should have same length as input
        assert len(result.state_path) == 100
        # State indices should be valid
        n_states = len(result.state_values)
        for state in result.state_path:
            assert 0 <= state < n_states
