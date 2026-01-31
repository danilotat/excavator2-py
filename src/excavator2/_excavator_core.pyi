"""
Type stubs for the EXCAVATOR2 C++ extension module.

This file provides type hints for IDE support and static type checking.
"""

from typing import List, overload

__version__: str

def hello() -> str:
    """Test function to verify C++ module works."""
    ...

def has_openmp() -> bool:
    """Check if OpenMP support is enabled."""
    ...

def get_num_threads() -> int:
    """Get number of OpenMP threads."""
    ...

def elnsum(x: float, y: float) -> float:
    """
    Log-sum-exp: compute log(exp(x) + exp(y)) in a numerically stable way.

    Args:
        x: First log-space value
        y: Second log-space value

    Returns:
        log(exp(x) + exp(y))
    """
    ...

def logsumexp(v: List[float]) -> float:
    """
    Log-sum-exp for a vector: compute log(sum(exp(v[i]))) in a numerically stable way.

    Args:
        v: Vector of log-space values

    Returns:
        log(sum(exp(v[i])))
    """
    ...

class hslm:
    """HSLM (Heterogeneous Shifting Level Model) segmentation algorithm module."""

    class HSLMParameters:
        """
        Parameters for the HSLM segmentation algorithm.

        Attributes:
            omega: Variance partitioning factor (0-1). Default: 0.1
            theta: Base transition probability. Default: 1e-5
            d_norm: Distance normalization factor. Default: 1e6
            step_eta: Distance step for eta computation (bp). Default: 200000
            n_states: Number of hidden states. Default: 21
            min_segment_size: Minimum segment size to keep. Default: 1
        """

        omega: float
        theta: float
        d_norm: float
        step_eta: float
        n_states: int
        min_segment_size: int

        def __init__(self) -> None: ...

    class HSLMResult:
        """
        Result of HSLM segmentation.

        Attributes:
            breakpoints: Indices where segments change (0-indexed)
            segment_means: Mean value for each segment
            state_path: State sequence from Viterbi algorithm
            state_values: The actual values corresponding to each state
            n_segments: Number of segments detected
            success: Whether the algorithm succeeded
            error_message: Error message if failed
        """

        breakpoints: List[int]
        segment_means: List[float]
        state_path: List[int]
        state_values: List[float]
        n_segments: int
        success: bool
        error_message: str

        def __init__(self) -> None: ...

    class HSLM:
        """
        HSLM (Heterogeneous Shifting Level Model) segmentation algorithm.

        This class implements an inhomogeneous HMM for detecting breakpoints
        in log2-ratio data from copy number analysis.

        Args:
            params: HSLMParameters object with algorithm settings
        """

        params: HSLMParameters

        def __init__(self, params: HSLMParameters = ...) -> None: ...
        def segment(self, log2_ratios: List[float], positions: List[int]) -> HSLMResult:
            """
            Run segmentation on a single chromosome.

            Args:
                log2_ratios: Vector of log2 ratio values (test/control)
                positions: Genomic positions corresponding to each log2 ratio

            Returns:
                HSLMResult containing breakpoints and segment information
            """
            ...

        def segment_multi(self, data_matrix: List[List[float]], positions: List[int]) -> HSLMResult:
            """
            Run segmentation on multiple sequences (multi-sample mode).

            Args:
                data_matrix: List of log2 ratio vectors, one per sample
                positions: Genomic positions (same for all samples)

            Returns:
                HSLMResult containing breakpoints
            """
            ...

    @staticmethod
    def segment(
        log2_ratios: List[float],
        positions: List[int],
        omega: float = 0.1,
        theta: float = 1e-5,
        step_eta: float = 200000.0,
        n_states: int = 21,
        min_segment_size: int = 1,
    ) -> HSLMResult:
        """
        Convenience function for single-call HSLM segmentation.

        Args:
            log2_ratios: Vector of log2 ratio values
            positions: Genomic positions
            omega: Variance partitioning factor (default: 0.1)
            theta: Base transition probability (default: 1e-5)
            step_eta: Distance step for eta (default: 200000)
            n_states: Number of hidden states (default: 21)
            min_segment_size: Minimum segment size to keep (default: 1)

        Returns:
            HSLMResult with segmentation results
        """
        ...

class fastcall:
    """FastCall CNV classification algorithm module."""

    class FastCallParameters:
        """
        Parameters for the FastCall CNV classification algorithm.

        Attributes:
            cellularity: Tumor purity (0.0-1.0). Default: 1.0
            thrd: Lower threshold for normal state. Default: 0.5
            thru: Upper threshold for normal state. Default: 0.35
            min_exons: Minimum exons per segment. Default: 4
            max_iterations: Maximum EM iterations. Default: 1000
            convergence: Convergence threshold. Default: 1e-5
        """

        cellularity: float
        thrd: float
        thru: float
        min_exons: int
        max_iterations: int
        convergence: float

        def __init__(self) -> None: ...

    class SegmentCall:
        """
        Copy number call for a single segment.

        Attributes:
            cn_call: Relative CN call (-2, -1, 0, +1, +2)
            absolute_cn: Absolute copy number (0, 1, 2, 3, 4+)
            probability: Posterior probability of the call
            state_index: Index of the most likely state (0-4)
            segment_mean: Mean log2 ratio of the segment
        """

        cn_call: int
        absolute_cn: int
        probability: float
        state_index: int
        segment_mean: float

        def __init__(self) -> None: ...

    class FastCallResult:
        """
        Result of FastCall CNV classification.

        Attributes:
            calls: Calls for each segment
            state_means: Fitted means for each state
            state_sds: Fitted standard deviations
            state_priors: Fitted prior probabilities
            iterations: Number of EM iterations
            converged: Whether EM converged
            success: Whether the algorithm succeeded
            error_message: Error message if failed
        """

        calls: List[SegmentCall]
        state_means: List[float]
        state_sds: List[float]
        state_priors: List[float]
        iterations: int
        converged: bool
        success: bool
        error_message: str

        def __init__(self) -> None: ...

    class FastCall:
        """
        FastCall CNV classification algorithm.

        Uses a 5-state Gaussian mixture model with EM to classify segments
        into copy number states.

        Args:
            params: FastCallParameters object with algorithm settings
        """

        params: FastCallParameters

        def __init__(self, params: FastCallParameters = ...) -> None: ...
        def call(self, segment_means: List[float], segment_sds: List[float] = []) -> FastCallResult:
            """
            Call copy number states for segments.

            Args:
                segment_means: Mean log2 ratios for each segment
                segment_sds: Standard deviations for each segment (optional)

            Returns:
                FastCallResult containing CN calls and probabilities
            """
            ...

    @staticmethod
    def call(
        segment_means: List[float],
        cellularity: float = 1.0,
        thrd: float = 0.5,
        thru: float = 0.35,
        max_iterations: int = 1000,
    ) -> FastCallResult:
        """
        Convenience function for single-call FastCall classification.

        Args:
            segment_means: Vector of segment mean log2 ratios
            cellularity: Tumor purity (default: 1.0)
            thrd: Lower threshold for normal state (default: 0.5)
            thru: Upper threshold for normal state (default: 0.35)
            max_iterations: Maximum EM iterations (default: 1000)

        Returns:
            FastCallResult with CN calls
        """
        ...
