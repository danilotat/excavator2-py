"""
Tests for the FastCall CNV classification algorithm.

These tests verify that the C++ implementation correctly ports the original
R implementation from EXCAVATOR2.
"""

import pytest
import numpy as np


class TestFastCallParameters:
    """Test FastCallParameters configuration."""

    def test_default_parameters(self):
        """Test default parameter values."""
        from excavator2._excavator_core.fastcall import FastCallParameters

        params = FastCallParameters()
        assert params.cellularity == 1.0
        assert params.thrd == 0.5
        assert params.thru == 0.35
        assert params.min_exons == 4
        assert params.max_iterations == 1000
        assert params.convergence == 1e-5

    def test_parameter_modification(self):
        """Test parameter values can be modified."""
        from excavator2._excavator_core.fastcall import FastCallParameters

        params = FastCallParameters()
        params.cellularity = 0.8
        params.thrd = 0.6
        params.max_iterations = 500

        assert params.cellularity == 0.8
        assert params.thrd == 0.6
        assert params.max_iterations == 500

    def test_parameter_repr(self):
        """Test parameter string representation."""
        from excavator2._excavator_core.fastcall import FastCallParameters

        params = FastCallParameters()
        repr_str = repr(params)
        assert "FastCallParameters" in repr_str
        assert "cellularity" in repr_str


class TestFastCallBasic:
    """Test basic FastCall functionality."""

    @pytest.fixture
    def mixed_segments(self):
        """Generate segments with mixed CN states."""
        np.random.seed(42)
        return np.concatenate([
            np.random.normal(0.0, 0.1, 50),    # Normal (CN=2)
            np.random.normal(-0.5, 0.1, 20),   # Deletion (CN=1)
            np.random.normal(0.58, 0.1, 15),   # Duplication (CN=3)
            np.random.normal(-1.5, 0.1, 10),   # Deep deletion (CN=0-1)
        ])

    @pytest.fixture
    def normal_segments(self):
        """Generate all normal segments."""
        np.random.seed(123)
        return np.random.normal(0.0, 0.1, 50)

    def test_call_returns_result(self):
        """Test that call() returns a FastCallResult object."""
        from excavator2._excavator_core.fastcall import FastCall, FastCallResult

        segment_means = [0.0] * 20
        caller = FastCall()
        result = caller.call(segment_means)

        assert isinstance(result, FastCallResult)
        assert result.success

    def test_call_mixed_segments(self, mixed_segments):
        """Test classification on mixed CN states."""
        from excavator2._excavator_core.fastcall import FastCall

        caller = FastCall()
        result = caller.call(mixed_segments.tolist())

        assert result.success
        assert result.converged
        assert len(result.calls) == len(mixed_segments)

        # Check that we detect multiple CN states
        cn_calls = [call.cn_call for call in result.calls]
        unique_calls = set(cn_calls)
        assert len(unique_calls) > 1, "Should detect multiple CN states"

    def test_call_detects_normal(self, normal_segments):
        """Test that normal segments are called as CN=2."""
        from excavator2._excavator_core.fastcall import FastCall

        caller = FastCall()
        result = caller.call(normal_segments.tolist())

        assert result.success

        # Most segments should be called as normal (CN=2, call=0)
        normal_count = sum(1 for call in result.calls if call.cn_call == 0)
        assert normal_count >= 40, f"Expected most segments to be normal, got {normal_count}/50"

    def test_call_detects_deletion(self):
        """Test that deletion segments are detected."""
        from excavator2._excavator_core.fastcall import FastCall

        np.random.seed(42)
        # Clear deletion signal
        segment_means = np.random.normal(-0.6, 0.1, 30)

        caller = FastCall()
        result = caller.call(segment_means.tolist())

        assert result.success

        # Should detect deletions (CN<2, call<0)
        deletion_count = sum(1 for call in result.calls if call.cn_call < 0)
        assert deletion_count >= 20, f"Expected deletions, got {deletion_count}/30"

    def test_call_detects_duplication(self):
        """Test that duplication segments are detected."""
        from excavator2._excavator_core.fastcall import FastCall

        np.random.seed(42)
        # Clear duplication signal (log2(3/2) ≈ 0.58)
        segment_means = np.random.normal(0.58, 0.1, 30)

        caller = FastCall()
        result = caller.call(segment_means.tolist())

        assert result.success

        # Should detect duplications (CN>2, call>0)
        duplication_count = sum(1 for call in result.calls if call.cn_call > 0)
        assert duplication_count >= 20, f"Expected duplications, got {duplication_count}/30"

    def test_empty_data_fails(self):
        """Test that empty data returns error."""
        from excavator2._excavator_core.fastcall import FastCall

        caller = FastCall()
        result = caller.call([])

        assert not result.success
        assert result.error_message != ""

    def test_single_segment(self):
        """Test classification of a single segment."""
        from excavator2._excavator_core.fastcall import FastCall

        caller = FastCall()
        result = caller.call([0.0])

        assert result.success
        assert len(result.calls) == 1
        # Should be called as normal
        assert result.calls[0].cn_call == 0


class TestSegmentCall:
    """Test SegmentCall structure."""

    def test_segment_call_attributes(self):
        """Test that SegmentCall has all required attributes."""
        from excavator2._excavator_core.fastcall import FastCall

        np.random.seed(42)
        segment_means = list(np.random.normal(0.0, 0.1, 10))

        caller = FastCall()
        result = caller.call(segment_means)

        call = result.calls[0]
        assert hasattr(call, 'cn_call')
        assert hasattr(call, 'absolute_cn')
        assert hasattr(call, 'probability')
        assert hasattr(call, 'state_index')
        assert hasattr(call, 'segment_mean')

    def test_segment_call_values(self):
        """Test that SegmentCall values are reasonable."""
        from excavator2._excavator_core.fastcall import FastCall

        np.random.seed(42)
        segment_means = list(np.random.normal(0.0, 0.1, 10))

        caller = FastCall()
        result = caller.call(segment_means)

        for call in result.calls:
            # CN call should be in {-2, -1, 0, +1, +2}
            assert call.cn_call in [-2, -1, 0, 1, 2]

            # Absolute CN should be in {0, 1, 2, 3, 4}
            assert call.absolute_cn in [0, 1, 2, 3, 4]

            # Probability should be in [0, 1]
            assert 0.0 <= call.probability <= 1.0

            # State index should be in [0, 4]
            assert 0 <= call.state_index <= 4

    def test_cn_call_consistency(self):
        """Test that cn_call and absolute_cn are consistent."""
        from excavator2._excavator_core.fastcall import FastCall

        np.random.seed(42)
        segment_means = list(np.random.normal(0.0, 0.3, 50))

        caller = FastCall()
        result = caller.call(segment_means)

        for call in result.calls:
            # cn_call = absolute_cn - 2
            expected_call = call.absolute_cn - 2
            assert call.cn_call == expected_call


class TestFastCallEM:
    """Test EM algorithm convergence and parameters."""

    def test_em_converges(self):
        """Test that EM algorithm converges."""
        from excavator2._excavator_core.fastcall import FastCall

        np.random.seed(42)
        # Generate data from mixture of gaussians
        segment_means = np.concatenate([
            np.random.normal(0.0, 0.1, 40),
            np.random.normal(-0.6, 0.1, 20),
            np.random.normal(0.6, 0.1, 20),
        ])

        caller = FastCall()
        result = caller.call(segment_means.tolist())

        assert result.success
        assert result.converged
        assert result.iterations > 0

    def test_fitted_parameters(self):
        """Test that fitted parameters are reasonable."""
        from excavator2._excavator_core.fastcall import FastCall

        np.random.seed(42)
        segment_means = list(np.random.normal(0.0, 0.1, 50))

        caller = FastCall()
        result = caller.call(segment_means)

        assert result.success

        # Should have 5 state means
        assert len(result.state_means) == 5
        assert len(result.state_sds) == 5
        assert len(result.state_priors) == 5

        # Means should be ordered (approximately)
        for i in range(4):
            assert result.state_means[i] < result.state_means[i + 1]

        # SDs should be positive
        for sd in result.state_sds:
            assert sd > 0

        # Priors should sum to ~1
        prior_sum = sum(result.state_priors)
        assert abs(prior_sum - 1.0) < 0.01

    def test_max_iterations_respected(self):
        """Test that max_iterations parameter is respected."""
        from excavator2._excavator_core.fastcall import FastCall, FastCallParameters

        np.random.seed(42)
        segment_means = list(np.random.normal(0.0, 0.1, 50))

        params = FastCallParameters()
        params.max_iterations = 2
        caller = FastCall(params)
        result = caller.call(segment_means)

        assert result.success
        assert result.iterations <= 2


class TestFastCallConvenienceFunction:
    """Test the convenience call() function."""

    def test_convenience_function(self):
        """Test the module-level call() function."""
        from excavator2._excavator_core import fastcall

        np.random.seed(42)
        segment_means = list(np.random.normal(0.0, 0.1, 30))

        result = fastcall.call(segment_means)

        assert result.success
        assert len(result.calls) == 30

    def test_convenience_function_with_params(self):
        """Test convenience function with custom parameters."""
        from excavator2._excavator_core import fastcall

        np.random.seed(42)
        segment_means = list(np.random.normal(0.0, 0.1, 30))

        result = fastcall.call(
            segment_means,
            cellularity=0.8,
            thrd=0.6,
            max_iterations=100
        )

        assert result.success


class TestFastCallProbabilities:
    """Test posterior probability calculations."""

    def test_probabilities_sum_to_one(self):
        """Test that probabilities are valid."""
        from excavator2._excavator_core.fastcall import FastCall

        np.random.seed(42)
        segment_means = list(np.random.normal(0.0, 0.1, 30))

        caller = FastCall()
        result = caller.call(segment_means)

        for call in result.calls:
            # Probability should be in [0, 1]
            assert 0.0 <= call.probability <= 1.0

    def test_high_probability_for_clear_signals(self):
        """Test that clear signals get high probability calls."""
        from excavator2._excavator_core.fastcall import FastCall

        # Very clear normal signal
        segment_means = [0.0] * 20

        caller = FastCall()
        result = caller.call(segment_means)

        # Most calls should have high probability
        high_prob_count = sum(1 for call in result.calls if call.probability > 0.9)
        assert high_prob_count >= 15, f"Expected high confidence calls, got {high_prob_count}/20"


class TestFastCallStateMapping:
    """Test CN state mapping."""

    def test_state_to_cn_mapping(self):
        """Test that states map correctly to CN values."""
        from excavator2._excavator_core.fastcall import FastCall

        np.random.seed(42)
        # Generate segments at each expected mean
        segment_means = [-3.0, -1.0, 0.0, 0.58, 1.0]

        caller = FastCall()
        result = caller.call(segment_means)

        # Check that states are assigned correctly
        # State 0 (mean=-3.0) → CN=0, call=-2
        # State 1 (mean=-1.0) → CN=1, call=-1
        # State 2 (mean=0.0) → CN=2, call=0
        # State 3 (mean=0.58) → CN=3, call=+1
        # State 4 (mean=1.0) → CN=4, call=+2

        for i, call in enumerate(result.calls):
            assert call.state_index == i, f"Segment {i} should be in state {i}"
            expected_cn = i
            expected_call = i - 2
            assert call.absolute_cn == expected_cn
            assert call.cn_call == expected_call


class TestFastCallEdgeCases:
    """Test edge cases and error handling."""

    def test_extreme_values(self):
        """Test handling of extreme log2 ratio values."""
        from excavator2._excavator_core.fastcall import FastCall

        segment_means = [-5.0, -2.0, 0.0, 2.0, 5.0]

        caller = FastCall()
        result = caller.call(segment_means)

        assert result.success
        assert len(result.calls) == 5

    def test_constant_values(self):
        """Test handling of constant segment means."""
        from excavator2._excavator_core.fastcall import FastCall

        segment_means = [0.0] * 50

        caller = FastCall()
        result = caller.call(segment_means)

        assert result.success
        # All should be called as normal
        normal_count = sum(1 for call in result.calls if call.cn_call == 0)
        assert normal_count == 50

    def test_single_outlier(self):
        """Test handling of a single outlier."""
        from excavator2._excavator_core.fastcall import FastCall

        np.random.seed(42)
        segment_means = list(np.random.normal(0.0, 0.1, 49)) + [10.0]

        caller = FastCall()
        result = caller.call(segment_means)

        assert result.success
        assert len(result.calls) == 50


class TestFastCallReproducibility:
    """Test reproducibility of results."""

    def test_deterministic_results(self):
        """Test that results are deterministic."""
        from excavator2._excavator_core.fastcall import FastCall

        np.random.seed(42)
        segment_means = list(np.random.normal(0.0, 0.2, 30))

        caller1 = FastCall()
        result1 = caller1.call(segment_means)

        caller2 = FastCall()
        result2 = caller2.call(segment_means)

        # Results should be identical
        assert result1.iterations == result2.iterations
        assert result1.converged == result2.converged

        for i in range(len(segment_means)):
            assert result1.calls[i].cn_call == result2.calls[i].cn_call
            assert abs(result1.calls[i].probability - result2.calls[i].probability) < 1e-10
