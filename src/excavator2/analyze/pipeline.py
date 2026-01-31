"""
CNV analysis pipeline.

This module provides the main analysis pipeline that combines log2 ratio
computation, HSLM segmentation, and FastCall CNV classification.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Dict, Union, Tuple
import logging

import numpy as np

from excavator2.analyze.segment import HSLMSegmenter, Segment, SegmentationResult, PreEstimatedParams
from excavator2.analyze.call import FastCallCaller, CNVCall, ClassificationResult, CopyNumberState
from excavator2.analyze.ratio import (
    Log2RatioResult,
    compute_log2_ratio,
    compute_log2_ratio_pooled,
    apply_cellularity_correction
)
from excavator2.prepare.normalize import NormalizationResult, load_normalized_counts


logger = logging.getLogger(__name__)


@dataclass
class CNVSegment:
    """A CNV segment with genomic coordinates and classification.

    Attributes:
        chrom: Chromosome name
        start: Start position (bp)
        end: End position (bp)
        start_idx: Start index in data array
        end_idx: End index in data array
        n_probes: Number of windows in segment
        segment_mean: Mean log2 ratio
        cn_call: Relative CN call (-2, -1, 0, +1, +2)
        absolute_cn: Absolute copy number (0, 1, 2, 3, 4)
        probability: Posterior probability of call
        state: CopyNumberState enum
        region_class: 'IN' or 'OUT' (majority)
    """
    chrom: str
    start: int
    end: int
    start_idx: int
    end_idx: int
    n_probes: int
    segment_mean: float
    cn_call: int
    absolute_cn: int
    probability: float
    state: CopyNumberState
    region_class: str = "IN"

    @property
    def length(self) -> int:
        """Segment length in base pairs."""
        return self.end - self.start

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
class ChromosomeResult:
    """Analysis results for a single chromosome.

    Attributes:
        chrom: Chromosome name
        segments: List of CNVSegments
        breakpoints: Raw breakpoint indices
        n_segments: Number of segments
        n_cnvs: Number of CNV calls
    """
    chrom: str
    segments: List[CNVSegment]
    breakpoints: List[int]
    n_segments: int
    n_cnvs: int


@dataclass
class WindowResult:
    """Per-window result data for output compatibility.

    Attributes:
        chrom: Chromosome name
        position: Window midpoint position
        start: Window start position
        end: Window end position
        log2_ratio: Raw log2 ratio for this window
        segment_mean: Mean log2 ratio of the segment this window belongs to
        region_class: 'IN' or 'OUT'
        segment_idx: Index of the segment this window belongs to
        cn_call: Relative CN call for the segment
        absolute_cn: Absolute copy number for the segment
        probability: Posterior probability of the call
    """
    chrom: str
    position: int
    start: int
    end: int
    log2_ratio: float
    segment_mean: float
    region_class: str
    segment_idx: int
    cn_call: int
    absolute_cn: int
    probability: float


@dataclass
class AnalysisResult:
    """Complete CNV analysis result for a sample pair.

    Attributes:
        test_sample: Test sample name
        control_sample: Control sample name (or 'PooledControl')
        segments: All CNV segments across all chromosomes
        chromosome_results: Results per chromosome
        n_segments: Total number of segments
        n_cnvs: Total number of CNV calls
        parameters: Algorithm parameters used
        window_results: Per-window results (for HSLM-style output)
    """
    test_sample: str
    control_sample: str
    segments: List[CNVSegment]
    chromosome_results: Dict[str, ChromosomeResult]
    n_segments: int
    n_cnvs: int
    parameters: dict = field(default_factory=dict)
    window_results: List[WindowResult] = field(default_factory=list)

    @property
    def n_deletions(self) -> int:
        """Number of deletion segments."""
        return sum(1 for s in self.segments if s.is_deletion)

    @property
    def n_gains(self) -> int:
        """Number of gain segments."""
        return sum(1 for s in self.segments if s.is_gain)

    def get_cnvs(self) -> List[CNVSegment]:
        """Get only CNV segments (excluding normal)."""
        return [s for s in self.segments if s.is_cnv]

    def get_deletions(self) -> List[CNVSegment]:
        """Get deletion segments only."""
        return [s for s in self.segments if s.is_deletion]

    def get_gains(self) -> List[CNVSegment]:
        """Get gain segments only."""
        return [s for s in self.segments if s.is_gain]


@dataclass
class AnalysisParameters:
    """Parameters for CNV analysis pipeline.

    HSLM parameters:
        omega: Variance partitioning factor (0-1). Default: 0.1
        theta: Base transition probability. Default: 1e-5
        step_eta: Distance step for eta computation (bp). Default: 200000
        n_states: Number of hidden states. Default: 21

    FastCall parameters:
        cellularity: Tumor purity (0-1). Default: 1.0
        thrd: Lower threshold for normal state. Default: 0.5
        thru: Upper threshold for normal state. Default: 0.35
        min_exons: Minimum exons per segment for calling. Default: 4
    """
    # HSLM parameters
    omega: float = 0.1
    theta: float = 1e-5
    step_eta: float = 200000.0
    n_states: int = 21

    # FastCall parameters
    cellularity: float = 1.0
    thrd: float = 0.5
    thru: float = 0.35
    min_exons: int = 4

    def to_dict(self) -> dict:
        """Convert parameters to dictionary."""
        return {
            'hslm': {
                'omega': self.omega,
                'theta': self.theta,
                'step_eta': self.step_eta,
                'n_states': self.n_states
            },
            'fastcall': {
                'cellularity': self.cellularity,
                'thrd': self.thrd,
                'thru': self.thru,
                'min_exons': self.min_exons
            }
        }


class CNVAnalyzer:
    """CNV analysis pipeline combining HSLM and FastCall.

    This class orchestrates the complete CNV detection workflow:
    1. Compute log2 ratios (test vs control)
    2. Run HSLM segmentation per chromosome
    3. Run FastCall classification on segments
    4. Combine and filter results

    Example:
        >>> analyzer = CNVAnalyzer(params)
        >>> result = analyzer.analyze_paired(test_data, control_data)
        >>> for segment in result.get_cnvs():
        ...     print(f"{segment.chrom}:{segment.start}-{segment.end} CN={segment.absolute_cn}")

    Args:
        params: AnalysisParameters object with algorithm settings
    """

    def __init__(self, params: Optional[AnalysisParameters] = None):
        if params is None:
            params = AnalysisParameters()
        self.params = params

        # Initialize algorithm components
        self._segmenter = HSLMSegmenter(
            omega=params.omega,
            theta=params.theta,
            step_eta=params.step_eta,
            n_states=params.n_states,
            min_segment_size=params.min_exons
        )

        self._caller = FastCallCaller(
            cellularity=params.cellularity,
            thrd=params.thrd,
            thru=params.thru,
            min_exons=params.min_exons
        )

    def analyze_paired(
        self,
        test_data: NormalizationResult,
        control_data: NormalizationResult,
        chromosomes: Optional[List[str]] = None
    ) -> AnalysisResult:
        """Analyze a paired test/control sample.

        Args:
            test_data: Normalized counts for test sample
            control_data: Normalized counts for control sample
            chromosomes: Chromosomes to analyze (default: all)

        Returns:
            AnalysisResult with CNV segments
        """
        logger.info(f"Starting paired analysis: {test_data.sample_name} vs {control_data.sample_name}")

        # Compute log2 ratios
        ratio_result = compute_log2_ratio(test_data, control_data)

        return self._analyze_ratios(ratio_result, chromosomes)

    def analyze_pooled(
        self,
        test_data: NormalizationResult,
        control_samples: List[NormalizationResult],
        chromosomes: Optional[List[str]] = None
    ) -> AnalysisResult:
        """Analyze a test sample against pooled controls.

        Args:
            test_data: Normalized counts for test sample
            control_samples: List of normalized counts for control samples
            chromosomes: Chromosomes to analyze (default: all)

        Returns:
            AnalysisResult with CNV segments
        """
        logger.info(f"Starting pooled analysis: {test_data.sample_name} vs {len(control_samples)} controls")

        # Compute log2 ratios against pooled control
        ratio_result = compute_log2_ratio_pooled(test_data, control_samples)

        return self._analyze_ratios(ratio_result, chromosomes)

    def _analyze_ratios(
        self,
        ratio_result: Log2RatioResult,
        chromosomes: Optional[List[str]] = None
    ) -> AnalysisResult:
        """Run segmentation and calling on log2 ratios.

        Args:
            ratio_result: Computed log2 ratios
            chromosomes: Chromosomes to analyze (default: all)

        Returns:
            AnalysisResult with CNV segments
        """
        if chromosomes is None:
            chromosomes = ratio_result.chromosomes

        all_segments = []
        all_window_results = []
        chromosome_results = {}

        # CRITICAL: Estimate HSLM parameters globally from ALL chromosomes' data
        # This matches the original R implementation where ParamEstSeq is called once
        # with all log2 ratios combined before per-chromosome segmentation.
        logger.info("Estimating HSLM parameters from global data...")
        global_params = self._segmenter.estimate_params(ratio_result.log2_ratios)
        if global_params.valid:
            logger.info(f"  Global parameters: smu={global_params.smu[0]:.4f}, sepsilon={global_params.sepsilon[0]:.4f}")
        else:
            logger.warning("  Failed to estimate global parameters, will use per-chromosome estimation")
            global_params = None

        for chrom in chromosomes:
            logger.info(f"Analyzing chromosome: {chrom}")

            # Get data for this chromosome
            chrom_data = ratio_result.get_chromosome_data(chrom)

            if chrom_data['n_windows'] == 0:
                logger.warning(f"  No windows for chromosome {chrom}, skipping")
                continue

            # Run analysis on this chromosome with global parameters
            chrom_result, window_results = self._analyze_chromosome(
                chrom=chrom,
                log2_ratios=chrom_data['log2_ratios'],
                positions=chrom_data['positions'],
                windows=chrom_data['windows'],
                in_target_mask=chrom_data['in_target_mask'],
                global_params=global_params
            )

            chromosome_results[chrom] = chrom_result
            all_segments.extend(chrom_result.segments)
            all_window_results.extend(window_results)

        # Summary statistics
        n_segments = len(all_segments)
        n_cnvs = sum(1 for s in all_segments if s.is_cnv)

        logger.info(f"Analysis complete: {n_segments} segments, {n_cnvs} CNVs")

        return AnalysisResult(
            test_sample=ratio_result.sample_name,
            control_sample=ratio_result.control_name,
            segments=all_segments,
            chromosome_results=chromosome_results,
            n_segments=n_segments,
            n_cnvs=n_cnvs,
            parameters=self.params.to_dict(),
            window_results=all_window_results
        )

    def _analyze_chromosome(
        self,
        chrom: str,
        log2_ratios: np.ndarray,
        positions: np.ndarray,
        windows: list,
        in_target_mask: np.ndarray,
        global_params: Optional[PreEstimatedParams] = None
    ) -> Tuple[ChromosomeResult, List[WindowResult]]:
        """Analyze a single chromosome.

        Args:
            chrom: Chromosome name
            log2_ratios: Log2 ratio values
            positions: Genomic positions
            windows: Window metadata
            in_target_mask: Boolean mask for IN-target windows
            global_params: Pre-estimated HSLM parameters (optional)

        Returns:
            Tuple of (ChromosomeResult, list of WindowResult)
        """
        # Store original log2 ratios for per-window output
        original_log2_ratios = log2_ratios.copy()

        # Apply cellularity correction if needed
        if self.params.cellularity < 1.0:
            log2_ratios = apply_cellularity_correction(
                log2_ratios, self.params.cellularity
            )

        # Step 1: HSLM Segmentation
        logger.info(f"  Running HSLM segmentation ({len(log2_ratios)} windows)")

        # Use global parameters if available (matching R behavior)
        if global_params is not None and global_params.valid:
            seg_result = self._segmenter.segment_with_params(log2_ratios, positions, global_params)
        else:
            seg_result = self._segmenter.segment(log2_ratios, positions)

        if not seg_result.success:
            logger.error(f"  Segmentation failed: {seg_result.error_message}")
            # Return empty results but still create window results with NA-like values
            window_results = self._create_window_results_no_seg(
                chrom, original_log2_ratios, windows, in_target_mask
            )
            return ChromosomeResult(
                chrom=chrom,
                segments=[],
                breakpoints=[],
                n_segments=0,
                n_cnvs=0
            ), window_results

        logger.info(f"  Found {seg_result.n_segments} segments")

        # Step 2: FastCall Classification
        if seg_result.n_segments > 0:
            segment_means = [s.mean for s in seg_result.segments]
            logger.info(f"  Running FastCall classification")
            call_result = self._caller.call(segment_means)

            if not call_result.success:
                logger.error(f"  Classification failed: {call_result.error_message}")
                # Return segments without calls
                chrom_result = self._segments_without_calls(chrom, seg_result, windows, in_target_mask)
                window_results = self._create_window_results(
                    chrom, original_log2_ratios, seg_result, None, windows, in_target_mask
                )
                return chrom_result, window_results

            logger.info(f"  EM converged in {call_result.iterations} iterations")
        else:
            call_result = None

        # Step 3: Combine segmentation and calling results
        segments = self._combine_results(
            chrom=chrom,
            seg_result=seg_result,
            call_result=call_result,
            windows=windows,
            positions=positions,
            in_target_mask=in_target_mask
        )

        # Filter by min_exons
        segments = self._filter_segments(segments)

        n_cnvs = sum(1 for s in segments if s.is_cnv)
        logger.info(f"  Final: {len(segments)} segments, {n_cnvs} CNVs")

        # Create per-window results
        window_results = self._create_window_results(
            chrom, original_log2_ratios, seg_result, call_result, windows, in_target_mask
        )

        return ChromosomeResult(
            chrom=chrom,
            segments=segments,
            breakpoints=seg_result.breakpoints,
            n_segments=len(segments),
            n_cnvs=n_cnvs
        ), window_results

    def _combine_results(
        self,
        chrom: str,
        seg_result: SegmentationResult,
        call_result: Optional[ClassificationResult],
        windows: list,
        positions: np.ndarray,
        in_target_mask: np.ndarray
    ) -> List[CNVSegment]:
        """Combine segmentation and classification results.

        Args:
            chrom: Chromosome name
            seg_result: HSLM segmentation result
            call_result: FastCall classification result
            windows: Window metadata
            positions: Genomic positions
            in_target_mask: Boolean mask for IN-target windows

        Returns:
            List of CNVSegment objects
        """
        segments = []

        for i, seg in enumerate(seg_result.segments):
            # Get CNV call if available
            if call_result and i < len(call_result.calls):
                call = call_result.calls[i]
                cn_call = call.cn_call
                absolute_cn = call.absolute_cn
                probability = call.probability
                state = call.state
            else:
                # Default to normal if no call available
                cn_call = 0
                absolute_cn = 2
                probability = 0.0
                state = CopyNumberState.NORMAL

            # Determine majority region class
            seg_in_mask = in_target_mask[seg.start_idx:seg.end_idx]
            region_class = 'IN' if np.sum(seg_in_mask) > len(seg_in_mask) / 2 else 'OUT'

            # Get genomic coordinates from windows
            seg_windows = windows[seg.start_idx:seg.end_idx]
            start_pos = seg_windows[0].start if seg_windows else seg.start_pos
            end_pos = seg_windows[-1].end if seg_windows else seg.end_pos

            segments.append(CNVSegment(
                chrom=chrom,
                start=start_pos,
                end=end_pos,
                start_idx=seg.start_idx,
                end_idx=seg.end_idx,
                n_probes=seg.n_probes,
                segment_mean=seg.mean,
                cn_call=cn_call,
                absolute_cn=absolute_cn,
                probability=probability,
                state=state,
                region_class=region_class
            ))

        return segments

    def _segments_without_calls(
        self,
        chrom: str,
        seg_result: SegmentationResult,
        windows: list,
        in_target_mask: np.ndarray
    ) -> ChromosomeResult:
        """Create result when classification fails."""
        segments = []

        for seg in seg_result.segments:
            seg_in_mask = in_target_mask[seg.start_idx:seg.end_idx]
            region_class = 'IN' if np.sum(seg_in_mask) > len(seg_in_mask) / 2 else 'OUT'

            seg_windows = windows[seg.start_idx:seg.end_idx]
            start_pos = seg_windows[0].start if seg_windows else seg.start_pos
            end_pos = seg_windows[-1].end if seg_windows else seg.end_pos

            segments.append(CNVSegment(
                chrom=chrom,
                start=start_pos,
                end=end_pos,
                start_idx=seg.start_idx,
                end_idx=seg.end_idx,
                n_probes=seg.n_probes,
                segment_mean=seg.mean,
                cn_call=0,
                absolute_cn=2,
                probability=0.0,
                state=CopyNumberState.NORMAL,
                region_class=region_class
            ))

        return ChromosomeResult(
            chrom=chrom,
            segments=segments,
            breakpoints=seg_result.breakpoints,
            n_segments=len(segments),
            n_cnvs=0
        )

    def _filter_segments(self, segments: List[CNVSegment]) -> List[CNVSegment]:
        """Filter segments by minimum probe count.

        Segments with fewer probes than min_exons are merged into
        adjacent segments or marked as normal.

        Args:
            segments: List of CNV segments

        Returns:
            Filtered list of segments
        """
        min_probes = self.params.min_exons

        # For now, just mark small segments as normal
        # More sophisticated merging could be added later
        filtered = []
        for seg in segments:
            if seg.n_probes < min_probes and seg.is_cnv:
                # Mark as normal
                filtered.append(CNVSegment(
                    chrom=seg.chrom,
                    start=seg.start,
                    end=seg.end,
                    start_idx=seg.start_idx,
                    end_idx=seg.end_idx,
                    n_probes=seg.n_probes,
                    segment_mean=seg.segment_mean,
                    cn_call=0,
                    absolute_cn=2,
                    probability=0.0,
                    state=CopyNumberState.NORMAL,
                    region_class=seg.region_class
                ))
            else:
                filtered.append(seg)

        return filtered

    def _create_window_results(
        self,
        chrom: str,
        log2_ratios: np.ndarray,
        seg_result: SegmentationResult,
        call_result: Optional[ClassificationResult],
        windows: list,
        in_target_mask: np.ndarray
    ) -> List[WindowResult]:
        """Create per-window results for HSLM-style output.

        This creates one WindowResult per input window, with:
        - log2_ratio: The original per-window log2 ratio
        - segment_mean: The mean of the segment this window belongs to
        - cn_call, absolute_cn, probability: From FastCall for the segment

        Args:
            chrom: Chromosome name
            log2_ratios: Original per-window log2 ratios
            seg_result: HSLM segmentation result
            call_result: FastCall classification result (may be None)
            windows: Window metadata
            in_target_mask: Boolean mask for IN-target windows

        Returns:
            List of WindowResult, one per window
        """
        window_results = []
        n_windows = len(log2_ratios)

        # Build a mapping from window index to segment info
        window_to_segment = np.zeros(n_windows, dtype=int)
        segment_means = np.zeros(n_windows)
        cn_calls = np.zeros(n_windows, dtype=int)
        absolute_cns = np.full(n_windows, 2, dtype=int)  # Default to normal
        probabilities = np.zeros(n_windows)

        for seg_idx, seg in enumerate(seg_result.segments):
            start_idx = seg.start_idx
            end_idx = seg.end_idx

            # Segment mean (median of log2 ratios in segment)
            seg_mean = seg.mean

            # Get call info if available
            if call_result and seg_idx < len(call_result.calls):
                call = call_result.calls[seg_idx]
                cn_call = call.cn_call
                absolute_cn = call.absolute_cn
                prob = call.probability
            else:
                cn_call = 0
                absolute_cn = 2
                prob = 0.0

            # Assign to all windows in this segment
            window_to_segment[start_idx:end_idx] = seg_idx
            segment_means[start_idx:end_idx] = seg_mean
            cn_calls[start_idx:end_idx] = cn_call
            absolute_cns[start_idx:end_idx] = absolute_cn
            probabilities[start_idx:end_idx] = prob

        # Create WindowResult for each window
        for i in range(n_windows):
            w = windows[i]
            region_class = 'IN' if in_target_mask[i] else 'OUT'

            window_results.append(WindowResult(
                chrom=chrom,
                position=w.position,
                start=w.start,
                end=w.end,
                log2_ratio=float(log2_ratios[i]),
                segment_mean=float(segment_means[i]),
                region_class=region_class,
                segment_idx=int(window_to_segment[i]),
                cn_call=int(cn_calls[i]),
                absolute_cn=int(absolute_cns[i]),
                probability=float(probabilities[i])
            ))

        return window_results

    def _create_window_results_no_seg(
        self,
        chrom: str,
        log2_ratios: np.ndarray,
        windows: list,
        in_target_mask: np.ndarray
    ) -> List[WindowResult]:
        """Create per-window results when segmentation failed.

        Uses NaN for segment_mean and default values for calls.

        Args:
            chrom: Chromosome name
            log2_ratios: Original per-window log2 ratios
            windows: Window metadata
            in_target_mask: Boolean mask for IN-target windows

        Returns:
            List of WindowResult, one per window
        """
        window_results = []

        for i, w in enumerate(windows):
            region_class = 'IN' if in_target_mask[i] else 'OUT'

            window_results.append(WindowResult(
                chrom=chrom,
                position=w.position,
                start=w.start,
                end=w.end,
                log2_ratio=float(log2_ratios[i]),
                segment_mean=float('nan'),  # No segmentation
                region_class=region_class,
                segment_idx=-1,  # No segment
                cn_call=0,
                absolute_cn=2,
                probability=0.0
            ))

        return window_results


def analyze_sample_pair(
    test_path: Union[str, Path],
    control_path: Union[str, Path],
    params: Optional[AnalysisParameters] = None,
    chromosomes: Optional[List[str]] = None
) -> AnalysisResult:
    """Convenience function to analyze a test/control pair.

    Args:
        test_path: Path to test sample NRC HDF5 file
        control_path: Path to control sample NRC HDF5 file
        params: Analysis parameters (default: AnalysisParameters())
        chromosomes: Chromosomes to analyze (default: all)

    Returns:
        AnalysisResult with CNV segments
    """
    test_data = load_normalized_counts(test_path)
    control_data = load_normalized_counts(control_path)

    analyzer = CNVAnalyzer(params)
    return analyzer.analyze_paired(test_data, control_data, chromosomes)
