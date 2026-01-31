"""
FastCall CNV classification algorithm wrapper.

This module provides a high-level Python interface to the C++ FastCall
CNV classification algorithm using Expectation-Maximization.
"""

from dataclasses import dataclass
from enum import IntEnum
from typing import List, Optional, Union
import numpy as np

from excavator2._excavator_core import fastcall


class CopyNumberState(IntEnum):
    """Copy number state enumeration.

    Values represent the relative copy number compared to diploid (2):
    - HOMOZYGOUS_DELETION: CN=0 (both copies lost)
    - HETEROZYGOUS_DELETION: CN=1 (one copy lost)
    - NORMAL: CN=2 (diploid, normal)
    - SINGLE_COPY_GAIN: CN=3 (one copy gained)
    - AMPLIFICATION: CN=4+ (high-level amplification)
    """

    HOMOZYGOUS_DELETION = -2
    HETEROZYGOUS_DELETION = -1
    NORMAL = 0
    SINGLE_COPY_GAIN = 1
    AMPLIFICATION = 2

    @classmethod
    def from_absolute_cn(cls, cn: int) -> "CopyNumberState":
        """Convert absolute copy number to CopyNumberState."""
        return cls(cn - 2)

    @property
    def absolute_cn(self) -> int:
        """Get absolute copy number (0, 1, 2, 3, 4)."""
        return self.value + 2

    @property
    def label(self) -> str:
        """Get human-readable label for this state."""
        labels = {
            -2: "Homozygous Deletion (CN=0)",
            -1: "Heterozygous Deletion (CN=1)",
            0: "Normal (CN=2)",
            1: "Single Copy Gain (CN=3)",
            2: "Amplification (CN=4+)",
        }
        return labels[self.value]


@dataclass
class CNVCall:
    """Copy number call for a single segment.

    Attributes:
        cn_call: Relative CN call (-2, -1, 0, +1, +2)
        absolute_cn: Absolute copy number (0, 1, 2, 3, 4)
        state: CopyNumberState enum value
        probability: Posterior probability of the call
        segment_mean: Mean log2 ratio of the segment
    """

    cn_call: int
    absolute_cn: int
    state: CopyNumberState
    probability: float
    segment_mean: float

    @property
    def is_cnv(self) -> bool:
        """Returns True if this is a CNV (not normal diploid)."""
        return self.cn_call != 0

    @property
    def is_deletion(self) -> bool:
        """Returns True if this is a deletion (CN < 2)."""
        return self.cn_call < 0

    @property
    def is_gain(self) -> bool:
        """Returns True if this is a gain (CN > 2)."""
        return self.cn_call > 0


@dataclass
class ClassificationResult:
    """Result of FastCall CNV classification.

    Attributes:
        calls: List of CNV calls, one per input segment
        state_means: Fitted mean log2 ratio for each CN state (5 values)
        state_sds: Fitted standard deviation for each CN state (5 values)
        state_priors: Fitted prior probability for each CN state (5 values)
        iterations: Number of EM iterations performed
        converged: Whether EM algorithm converged
        success: Whether classification succeeded
        error_message: Error message if failed
    """

    calls: List[CNVCall]
    state_means: List[float]
    state_sds: List[float]
    state_priors: List[float]
    iterations: int
    converged: bool
    success: bool
    error_message: str = ""

    @property
    def n_segments(self) -> int:
        """Number of segments classified."""
        return len(self.calls)

    @property
    def n_cnvs(self) -> int:
        """Number of segments with CNV (non-normal)."""
        return sum(1 for c in self.calls if c.is_cnv)

    @property
    def n_deletions(self) -> int:
        """Number of deletion segments."""
        return sum(1 for c in self.calls if c.is_deletion)

    @property
    def n_gains(self) -> int:
        """Number of gain segments."""
        return sum(1 for c in self.calls if c.is_gain)

    def get_cnvs(self) -> List[CNVCall]:
        """Get only CNV calls (excluding normal segments)."""
        return [c for c in self.calls if c.is_cnv]

    def get_deletions(self) -> List[CNVCall]:
        """Get only deletion calls."""
        return [c for c in self.calls if c.is_deletion]

    def get_gains(self) -> List[CNVCall]:
        """Get only gain calls."""
        return [c for c in self.calls if c.is_gain]


class FastCallCaller:
    """High-level wrapper for FastCall C++ implementation.

    FastCall uses a 5-state Gaussian mixture model with Expectation-
    Maximization to classify segments into copy number states:
    - CN=0: Homozygous deletion (mean ~ -3.0)
    - CN=1: Heterozygous deletion (mean ~ -1.0)
    - CN=2: Normal/diploid (mean ~ 0.0)
    - CN=3: Single copy gain (mean ~ 0.58)
    - CN=4+: Amplification (mean ~ 1.0)

    Example:
        >>> caller = FastCallCaller(cellularity=0.8)
        >>> result = caller.call(segment_means)
        >>> for i, call in enumerate(result.calls):
        ...     if call.is_cnv:
        ...         print(f"Segment {i}: {call.state.label}, prob={call.probability:.3f}")

    Args:
        cellularity: Tumor purity (0.0-1.0). Lower values adjust expected
            log2 ratios for admixture with normal cells. Default: 1.0
        thrd: Lower threshold for normal state (deletion boundary).
            Default: 0.5
        thru: Upper threshold for normal state (duplication boundary).
            Default: 0.35
        min_exons: Minimum exons per segment for calling. Default: 4
        max_iterations: Maximum EM iterations. Default: 1000
        convergence: Convergence threshold for EM. Default: 1e-5
    """

    def __init__(
        self,
        cellularity: float = 1.0,
        thrd: float = 0.5,
        thru: float = 0.35,
        min_exons: int = 4,
        max_iterations: int = 1000,
        convergence: float = 1e-5,
    ):
        if not 0 < cellularity <= 1:
            raise ValueError(f"cellularity must be in (0, 1], got {cellularity}")
        if thrd <= 0:
            raise ValueError(f"thrd must be positive, got {thrd}")
        if thru <= 0:
            raise ValueError(f"thru must be positive, got {thru}")

        self._params = fastcall.FastCallParameters()
        self._params.cellularity = cellularity
        self._params.thrd = thrd
        self._params.thru = thru
        self._params.min_exons = min_exons
        self._params.max_iterations = max_iterations
        self._params.convergence = convergence
        self._caller = fastcall.FastCall(self._params)

    @property
    def cellularity(self) -> float:
        """Tumor purity/cellularity."""
        return self._params.cellularity

    @cellularity.setter
    def cellularity(self, value: float):
        if not 0 < value <= 1:
            raise ValueError(f"cellularity must be in (0, 1], got {value}")
        self._params.cellularity = value
        self._caller = fastcall.FastCall(self._params)

    @property
    def max_iterations(self) -> int:
        """Maximum EM iterations."""
        return self._params.max_iterations

    @max_iterations.setter
    def max_iterations(self, value: int):
        if value <= 0:
            raise ValueError(f"max_iterations must be positive, got {value}")
        self._params.max_iterations = value
        self._caller = fastcall.FastCall(self._params)

    def call(
        self,
        segment_means: Union[np.ndarray, List[float]],
        segment_sds: Optional[Union[np.ndarray, List[float]]] = None,
    ) -> ClassificationResult:
        """Call copy number states for segments.

        Args:
            segment_means: Mean log2 ratios for each segment
            segment_sds: Standard deviations for each segment (optional)

        Returns:
            ClassificationResult with CN calls and EM fit statistics

        Raises:
            ValueError: If inputs are invalid
        """
        # Convert to numpy for validation
        means = np.asarray(segment_means, dtype=np.float64)

        if len(means) == 0:
            raise ValueError("segment_means cannot be empty")

        # Check for NaN/Inf
        if np.any(~np.isfinite(means)):
            raise ValueError("segment_means contains NaN or Inf values")

        # Convert segment_sds if provided
        sds_list = []
        if segment_sds is not None:
            sds = np.asarray(segment_sds, dtype=np.float64)
            if len(sds) != len(means):
                raise ValueError(
                    f"Length mismatch: segment_means ({len(means)}) vs " f"segment_sds ({len(sds)})"
                )
            sds_list = sds.tolist()

        # Call C++ implementation
        result = self._caller.call(means.tolist(), sds_list)

        if not result.success:
            return ClassificationResult(
                calls=[],
                state_means=[],
                state_sds=[],
                state_priors=[],
                iterations=0,
                converged=False,
                success=False,
                error_message=result.error_message,
            )

        # Convert C++ SegmentCall objects to Python CNVCall objects
        calls = []
        for c in result.calls:
            calls.append(
                CNVCall(
                    cn_call=c.cn_call,
                    absolute_cn=c.absolute_cn,
                    state=CopyNumberState(c.cn_call),
                    probability=c.probability,
                    segment_mean=c.segment_mean,
                )
            )

        return ClassificationResult(
            calls=calls,
            state_means=list(result.state_means),
            state_sds=list(result.state_sds),
            state_priors=list(result.state_priors),
            iterations=result.iterations,
            converged=result.converged,
            success=True,
        )


def call(
    segment_means: Union[np.ndarray, List[float]],
    cellularity: float = 1.0,
    thrd: float = 0.5,
    thru: float = 0.35,
    max_iterations: int = 1000,
) -> ClassificationResult:
    """Convenience function for single-call FastCall classification.

    This is equivalent to creating a FastCallCaller and calling call(),
    but more convenient for one-off usage.

    Args:
        segment_means: Mean log2 ratios for each segment
        cellularity: Tumor purity (default: 1.0)
        thrd: Lower threshold for normal state (default: 0.5)
        thru: Upper threshold for normal state (default: 0.35)
        max_iterations: Maximum EM iterations (default: 1000)

    Returns:
        ClassificationResult with CN calls

    Example:
        >>> result = call(segment_means, cellularity=0.8)
        >>> print(f"Found {result.n_cnvs} CNVs")
    """
    caller = FastCallCaller(
        cellularity=cellularity, thrd=thrd, thru=thru, max_iterations=max_iterations
    )
    return caller.call(segment_means)
