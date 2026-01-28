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

from excavator2.analyze.segment import HSLMSegmenter, Segment, SegmentationResult
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
    """
    test_sample: str
    control_sample: str
    segments: List[CNVSegment]
    chromosome_results: Dict[str, ChromosomeResult]
    n_segments: int
    n_cnvs: int
    parameters: dict = field(default_factory=dict)

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
            min_segment_size=1
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
        chromosome_results = {}

        for chrom in chromosomes:
            logger.info(f"Analyzing chromosome: {chrom}")

            # Get data for this chromosome
            chrom_data = ratio_result.get_chromosome_data(chrom)

            if chrom_data['n_windows'] == 0:
                logger.warning(f"  No windows for chromosome {chrom}, skipping")
                continue

            # Run analysis on this chromosome
            chrom_result = self._analyze_chromosome(
                chrom=chrom,
                log2_ratios=chrom_data['log2_ratios'],
                positions=chrom_data['positions'],
                windows=chrom_data['windows'],
                in_target_mask=chrom_data['in_target_mask']
            )

            chromosome_results[chrom] = chrom_result
            all_segments.extend(chrom_result.segments)

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
            parameters=self.params.to_dict()
        )

    def _analyze_chromosome(
        self,
        chrom: str,
        log2_ratios: np.ndarray,
        positions: np.ndarray,
        windows: list,
        in_target_mask: np.ndarray
    ) -> ChromosomeResult:
        """Analyze a single chromosome.

        Args:
            chrom: Chromosome name
            log2_ratios: Log2 ratio values
            positions: Genomic positions
            windows: Window metadata
            in_target_mask: Boolean mask for IN-target windows

        Returns:
            ChromosomeResult with segments
        """
        # Apply cellularity correction if needed
        if self.params.cellularity < 1.0:
            log2_ratios = apply_cellularity_correction(
                log2_ratios, self.params.cellularity
            )

        # Step 1: HSLM Segmentation
        logger.info(f"  Running HSLM segmentation ({len(log2_ratios)} windows)")
        seg_result = self._segmenter.segment(log2_ratios, positions)

        if not seg_result.success:
            logger.error(f"  Segmentation failed: {seg_result.error_message}")
            return ChromosomeResult(
                chrom=chrom,
                segments=[],
                breakpoints=[],
                n_segments=0,
                n_cnvs=0
            )

        logger.info(f"  Found {seg_result.n_segments} segments")

        # Step 2: FastCall Classification
        if seg_result.n_segments > 0:
            segment_means = [s.mean for s in seg_result.segments]
            logger.info(f"  Running FastCall classification")
            call_result = self._caller.call(segment_means)

            if not call_result.success:
                logger.error(f"  Classification failed: {call_result.error_message}")
                # Return segments without calls
                return self._segments_without_calls(chrom, seg_result, windows, in_target_mask)

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

        return ChromosomeResult(
            chrom=chrom,
            segments=segments,
            breakpoints=seg_result.breakpoints,
            n_segments=len(segments),
            n_cnvs=n_cnvs
        )

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
