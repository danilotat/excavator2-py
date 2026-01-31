"""
HSLM segmentation algorithm wrapper.

This module provides a high-level Python interface to the C++ HSLM
(Heterogeneous Shifting Level Model) segmentation algorithm.
"""

from dataclasses import dataclass
from typing import List, Optional, Dict, Any, Union
import numpy as np

from excavator2._excavator_core import hslm


@dataclass
class Segment:
    """A genomic segment detected by HSLM.

    Attributes:
        start_idx: Start index in the original data (0-indexed, inclusive)
        end_idx: End index in the original data (0-indexed, exclusive)
        start_pos: Genomic start position (bp)
        end_pos: Genomic end position (bp)
        mean: Mean log2 ratio of the segment
        n_probes: Number of probes/windows in the segment
    """

    start_idx: int
    end_idx: int
    start_pos: int
    end_pos: int
    mean: float
    n_probes: int


@dataclass
class PreEstimatedParams:
    """Pre-estimated parameters for HSLM segmentation.

    These can be computed once from all chromosomes' data,
    then used for per-chromosome segmentation (matching R behavior).

    Attributes:
        mi: Mean for each sequence (typically 0)
        smu: State standard deviation
        sepsilon: Noise standard deviation
        muk: State means matrix (n_sequences x n_states)
        valid: Whether parameters are valid
    """

    mi: List[float]
    smu: List[float]
    sepsilon: List[float]
    muk: List[List[float]]
    valid: bool = True


@dataclass
class SegmentationResult:
    """Result of HSLM segmentation.

    Attributes:
        segments: List of detected segments
        breakpoints: Raw breakpoint indices from HSLM
        state_path: HMM state path from Viterbi algorithm
        n_segments: Number of segments detected
        success: Whether segmentation succeeded
        error_message: Error message if failed
    """

    segments: List[Segment]
    breakpoints: List[int]
    state_path: List[int]
    n_segments: int
    success: bool
    error_message: str = ""


class HSLMSegmenter:
    """High-level wrapper for HSLM C++ implementation.

    The HSLM (Heterogeneous Shifting Level Model) algorithm detects
    breakpoints in log2 ratio profiles using an inhomogeneous Hidden
    Markov Model with distance-dependent transition probabilities.

    Example:
        >>> segmenter = HSLMSegmenter(omega=0.1, theta=1e-5)
        >>> result = segmenter.segment(log2_ratios, positions)
        >>> for seg in result.segments:
        ...     print(f"Segment: {seg.start_pos}-{seg.end_pos}, mean={seg.mean:.3f}")

    Args:
        omega: Variance partitioning factor (0-1). Higher values give more
            weight to within-state variance. Default: 0.1
        theta: Base transition probability. Lower values favor fewer
            breakpoints. Default: 1e-5
        step_eta: Distance step for eta computation (bp). Controls how
            distance affects transition probability. Default: 200000
        n_states: Number of hidden states for discretization. Default: 21
        min_segment_size: Minimum segment size to keep. Default: 1
    """

    def __init__(
        self,
        omega: float = 0.1,
        theta: float = 1e-5,
        step_eta: float = 200000.0,
        n_states: int = 21,
        min_segment_size: int = 1,
    ):
        self._params = hslm.HSLMParameters()
        self._params.omega = omega
        self._params.theta = theta
        self._params.step_eta = step_eta
        self._params.n_states = n_states
        self._params.min_segment_size = min_segment_size
        self._hslm = hslm.HSLM(self._params)

    @property
    def omega(self) -> float:
        """Variance partitioning factor."""
        return self._params.omega

    @omega.setter
    def omega(self, value: float):
        if not 0 <= value <= 1:
            raise ValueError(f"omega must be in [0, 1], got {value}")
        self._params.omega = value
        self._hslm = hslm.HSLM(self._params)

    @property
    def theta(self) -> float:
        """Base transition probability."""
        return self._params.theta

    @theta.setter
    def theta(self, value: float):
        if value <= 0:
            raise ValueError(f"theta must be positive, got {value}")
        self._params.theta = value
        self._hslm = hslm.HSLM(self._params)

    @property
    def n_states(self) -> int:
        """Number of hidden states."""
        return self._params.n_states

    @property
    def min_segment_size(self) -> int:
        """Minimum segment size."""
        return self._params.min_segment_size

    def segment(
        self, log2_ratios: Union[np.ndarray, List[float]], positions: Union[np.ndarray, List[int]]
    ) -> SegmentationResult:
        """Run HSLM segmentation on a single chromosome/region.

        Args:
            log2_ratios: Array of log2 ratio values (test/control)
            positions: Genomic positions corresponding to each log2 ratio

        Returns:
            SegmentationResult with detected segments and metadata

        Raises:
            ValueError: If inputs have mismatched lengths or are empty
        """
        # Convert to numpy arrays for validation
        ratios = np.asarray(log2_ratios, dtype=np.float64)
        pos = np.asarray(positions, dtype=np.int64)

        # Input validation
        if len(ratios) == 0:
            raise ValueError("log2_ratios cannot be empty")
        if len(ratios) != len(pos):
            raise ValueError(
                f"Length mismatch: log2_ratios ({len(ratios)}) vs positions ({len(pos)})"
            )

        # Check for NaN/Inf values
        if np.any(~np.isfinite(ratios)):
            raise ValueError("log2_ratios contains NaN or Inf values")

        # Call C++ implementation
        result = self._hslm.segment(ratios.tolist(), pos.tolist())

        if not result.success:
            return SegmentationResult(
                segments=[],
                breakpoints=[],
                state_path=[],
                n_segments=0,
                success=False,
                error_message=result.error_message,
            )

        # Convert breakpoints to Segment objects
        segments = self._breakpoints_to_segments(
            result.breakpoints, result.segment_means, ratios, pos
        )

        return SegmentationResult(
            segments=segments,
            breakpoints=list(result.breakpoints),
            state_path=list(result.state_path),
            n_segments=result.n_segments,
            success=True,
        )

    def estimate_params(self, log2_ratios: Union[np.ndarray, List[float]]) -> PreEstimatedParams:
        """Estimate HSLM parameters from data without running segmentation.

        This allows computing parameters once from all chromosomes' data,
        then using them for per-chromosome segmentation (matching R behavior).

        Args:
            log2_ratios: Array of log2 ratio values (all chromosomes combined)

        Returns:
            PreEstimatedParams containing mi, smu, sepsilon, muk
        """
        ratios = np.asarray(log2_ratios, dtype=np.float64)

        if len(ratios) == 0:
            return PreEstimatedParams(mi=[], smu=[], sepsilon=[], muk=[], valid=False)

        # Call C++ implementation with single-sample matrix
        data_matrix = [ratios.tolist()]
        result = self._hslm.estimate_params(data_matrix)

        return PreEstimatedParams(
            mi=list(result.mi),
            smu=list(result.smu),
            sepsilon=list(result.sepsilon),
            muk=[list(row) for row in result.muk],
            valid=result.valid,
        )

    def segment_with_params(
        self,
        log2_ratios: Union[np.ndarray, List[float]],
        positions: Union[np.ndarray, List[int]],
        params: PreEstimatedParams,
    ) -> SegmentationResult:
        """Run segmentation with pre-estimated parameters.

        This matches the original R behavior where parameters are estimated
        globally from all data, then used for per-chromosome segmentation.

        Args:
            log2_ratios: Array of log2 ratio values
            positions: Genomic positions
            params: Pre-estimated parameters from estimate_params()

        Returns:
            SegmentationResult with detected segments
        """
        ratios = np.asarray(log2_ratios, dtype=np.float64)
        pos = np.asarray(positions, dtype=np.int64)

        if len(ratios) == 0:
            raise ValueError("log2_ratios cannot be empty")
        if len(ratios) != len(pos):
            raise ValueError(
                f"Length mismatch: log2_ratios ({len(ratios)}) vs positions ({len(pos)})"
            )

        if np.any(~np.isfinite(ratios)):
            raise ValueError("log2_ratios contains NaN or Inf values")

        # Convert PreEstimatedParams to C++ format
        cpp_params = hslm.PreEstimatedParams()
        cpp_params.mi = params.mi
        cpp_params.smu = params.smu
        cpp_params.sepsilon = params.sepsilon
        cpp_params.muk = params.muk
        cpp_params.valid = params.valid

        result = self._hslm.segment_with_params(ratios.tolist(), pos.tolist(), cpp_params)

        if not result.success:
            return SegmentationResult(
                segments=[],
                breakpoints=[],
                state_path=[],
                n_segments=0,
                success=False,
                error_message=result.error_message,
            )

        segments = self._breakpoints_to_segments(
            result.breakpoints, result.segment_means, ratios, pos
        )

        return SegmentationResult(
            segments=segments,
            breakpoints=list(result.breakpoints),
            state_path=list(result.state_path),
            n_segments=result.n_segments,
            success=True,
        )

    def segment_multi(
        self,
        data_matrix: Union[np.ndarray, List[List[float]]],
        positions: Union[np.ndarray, List[int]],
    ) -> SegmentationResult:
        """Run segmentation on multiple samples (multi-sample mode).

        In multi-sample mode, breakpoints are detected jointly across
        all samples, which can improve detection for shared events.

        Args:
            data_matrix: 2D array where each row is a sample's log2 ratios
            positions: Genomic positions (same for all samples)

        Returns:
            SegmentationResult with detected segments

        Raises:
            ValueError: If inputs are invalid
        """
        # Convert to numpy array
        data = np.asarray(data_matrix, dtype=np.float64)
        pos = np.asarray(positions, dtype=np.int64)

        if data.ndim != 2:
            raise ValueError(f"data_matrix must be 2D, got {data.ndim}D")
        if data.shape[1] != len(pos):
            raise ValueError(f"Column count ({data.shape[1]}) != positions length ({len(pos)})")

        # Convert to list of lists for C++
        data_list = [row.tolist() for row in data]

        result = self._hslm.segment_multi(data_list, pos.tolist())

        if not result.success:
            return SegmentationResult(
                segments=[],
                breakpoints=[],
                state_path=[],
                n_segments=0,
                success=False,
                error_message=result.error_message,
            )

        # For multi-sample, use mean across samples for segment means
        mean_ratios = np.mean(data, axis=0)
        segments = self._breakpoints_to_segments(
            result.breakpoints, result.segment_means, mean_ratios, pos
        )

        return SegmentationResult(
            segments=segments,
            breakpoints=list(result.breakpoints),
            state_path=list(result.state_path),
            n_segments=result.n_segments,
            success=True,
        )

    def _breakpoints_to_segments(
        self,
        breakpoints: List[int],
        segment_means: List[float],
        ratios: np.ndarray,
        positions: np.ndarray,
    ) -> List[Segment]:
        """Convert breakpoint indices to Segment objects."""
        segments = []
        n = len(ratios)

        # Add endpoint to breakpoints if not present
        bp_list = list(breakpoints)
        if not bp_list or bp_list[-1] != n:
            bp_list.append(n)

        start_idx = 0
        for i, end_idx in enumerate(bp_list):
            if end_idx <= start_idx:
                continue

            # Get segment mean from HSLM result if available, else compute
            if i < len(segment_means):
                mean = segment_means[i]
            else:
                mean = float(np.mean(ratios[start_idx:end_idx]))

            segments.append(
                Segment(
                    start_idx=start_idx,
                    end_idx=end_idx,
                    start_pos=int(positions[start_idx]),
                    end_pos=int(positions[end_idx - 1]),
                    mean=mean,
                    n_probes=end_idx - start_idx,
                )
            )

            start_idx = end_idx

        return segments


def segment(
    log2_ratios: Union[np.ndarray, List[float]],
    positions: Union[np.ndarray, List[int]],
    omega: float = 0.1,
    theta: float = 1e-5,
    step_eta: float = 200000.0,
    n_states: int = 21,
    min_segment_size: int = 1,
) -> SegmentationResult:
    """Convenience function for single-call HSLM segmentation.

    This is equivalent to creating an HSLMSegmenter and calling segment(),
    but more convenient for one-off usage.

    Args:
        log2_ratios: Array of log2 ratio values
        positions: Genomic positions
        omega: Variance partitioning factor (default: 0.1)
        theta: Base transition probability (default: 1e-5)
        step_eta: Distance step for eta (default: 200000)
        n_states: Number of hidden states (default: 21)
        min_segment_size: Minimum segment size (default: 1)

    Returns:
        SegmentationResult with detected segments

    Example:
        >>> result = segment(log2_ratios, positions, omega=0.1)
        >>> print(f"Found {result.n_segments} segments")
    """
    segmenter = HSLMSegmenter(
        omega=omega,
        theta=theta,
        step_eta=step_eta,
        n_states=n_states,
        min_segment_size=min_segment_size,
    )
    return segmenter.segment(log2_ratios, positions)
