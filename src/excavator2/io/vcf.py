"""
VCF output writer for CNV calls.

This module provides functions for writing CNV calls to VCF format,
following the VCF 4.2 specification for structural variants.
"""

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Union, TextIO
import logging

import numpy as np

from excavator2.analyze.pipeline import CNVSegment, AnalysisResult
from excavator2.analyze.call import CopyNumberState


logger = logging.getLogger(__name__)


# VCF header template
VCF_HEADER_TEMPLATE = """##fileformat=VCFv4.2
##fileDate={date}
##source=EXCAVATOR2
##reference={reference}
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number estimate">
##INFO=<ID=CNREL,Number=1,Type=Integer,Description="Relative copy number (-2 to +2)">
##INFO=<ID=LOG2R,Number=1,Type=Float,Description="Log2 ratio of segment">
##INFO=<ID=NPROBES,Number=1,Type=Integer,Description="Number of probes/windows in segment">
##INFO=<ID=PROB,Number=1,Type=Float,Description="Posterior probability of call">
##INFO=<ID=CLASS,Number=1,Type=String,Description="Region class (IN or OUT target)">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=CNV,Description="Copy number variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number quality (probability)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"""


def _get_sv_type(segment: CNVSegment) -> str:
    """Get structural variant type for a segment."""
    if segment.cn_call < 0:
        return "DEL"
    elif segment.cn_call > 0:
        return "DUP"
    else:
        return "CNV"


def _get_alt_allele(segment: CNVSegment) -> str:
    """Get ALT allele representation for a segment."""
    if segment.cn_call < 0:
        return "<DEL>"
    elif segment.cn_call > 0:
        return "<DUP>"
    else:
        return "<CNV>"


def _segment_to_vcf_line(
    segment: CNVSegment,
    variant_id: str,
    sample_name: str
) -> str:
    """Convert a CNV segment to a VCF line.

    Args:
        segment: CNVSegment to convert
        variant_id: Unique variant identifier
        sample_name: Sample name

    Returns:
        Tab-separated VCF line
    """
    chrom = segment.chrom
    pos = segment.start
    ref = "N"  # Reference base (placeholder)
    alt = _get_alt_allele(segment)
    qual = "."  # Quality score
    filt = "PASS" if segment.is_cnv else "."

    # INFO field
    sv_type = _get_sv_type(segment)
    sv_len = segment.end - segment.start
    if segment.is_deletion:
        sv_len = -sv_len  # Deletions have negative SVLEN

    info_parts = [
        f"SVTYPE={sv_type}",
        f"END={segment.end}",
        f"SVLEN={sv_len}",
        f"CN={segment.absolute_cn}",
        f"CNREL={segment.cn_call}",
        f"LOG2R={segment.segment_mean:.4f}",
        f"NPROBES={segment.n_probes}",
        f"PROB={segment.probability:.4f}",
        f"CLASS={segment.region_class}"
    ]
    info = ";".join(info_parts)

    # FORMAT and sample fields
    fmt = "GT:CN:CNQ"

    # Genotype based on copy number
    if segment.absolute_cn == 0:
        gt = "1/1"  # Homozygous deletion
    elif segment.absolute_cn == 1:
        gt = "0/1"  # Heterozygous deletion
    elif segment.absolute_cn == 2:
        gt = "0/0"  # Normal
    elif segment.absolute_cn == 3:
        gt = "0/1"  # Single copy gain
    else:
        gt = "1/1"  # Amplification

    sample_data = f"{gt}:{segment.absolute_cn}:{segment.probability:.4f}"

    return f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\t{fmt}\t{sample_data}"


def write_vcf(
    result: AnalysisResult,
    output_path: Union[str, Path],
    reference: str = "unknown",
    cnv_only: bool = True
) -> None:
    """Write CNV calls to VCF file.

    Args:
        result: AnalysisResult containing CNV segments
        output_path: Path to output VCF file
        reference: Reference genome name (for header)
        cnv_only: If True, only write CNV calls (exclude normal segments)
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Writing VCF to {output_path}")

    with open(output_path, 'w') as f:
        # Write header
        header = VCF_HEADER_TEMPLATE.format(
            date=datetime.now().strftime("%Y%m%d"),
            reference=reference,
            sample=result.test_sample
        )
        f.write(header + "\n")

        # Write variants
        variant_count = 0
        for i, segment in enumerate(result.segments):
            if cnv_only and not segment.is_cnv:
                continue

            variant_id = f"EXCAVATOR2_{i+1}"
            line = _segment_to_vcf_line(segment, variant_id, result.test_sample)
            f.write(line + "\n")
            variant_count += 1

    logger.info(f"Wrote {variant_count} variants to VCF")


def write_vcf_regions(
    result: AnalysisResult,
    output_path: Union[str, Path],
    reference: str = "unknown"
) -> None:
    """Write region-level CNV calls to VCF (CNVs only).

    This writes only CNV calls (not normal segments), suitable for
    downstream analysis of detected variants.

    Args:
        result: AnalysisResult containing CNV segments
        output_path: Path to output VCF file
        reference: Reference genome name
    """
    write_vcf(result, output_path, reference, cnv_only=True)


def write_vcf_windows(
    result: AnalysisResult,
    output_path: Union[str, Path],
    reference: str = "unknown"
) -> None:
    """Write window-level calls to VCF (all segments).

    This writes all segments including normal regions, useful for
    visualization and QC purposes.

    Args:
        result: AnalysisResult containing CNV segments
        output_path: Path to output VCF file
        reference: Reference genome name
    """
    write_vcf(result, output_path, reference, cnv_only=False)


@dataclass
class BEDRecord:
    """BED format record for CNV output."""
    chrom: str
    start: int
    end: int
    name: str
    score: float
    strand: str = "."

    def to_line(self) -> str:
        """Convert to BED line."""
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.name}\t{self.score}\t{self.strand}"


def write_bed(
    result: AnalysisResult,
    output_path: Union[str, Path],
    cnv_only: bool = True
) -> None:
    """Write CNV calls to BED format.

    The BED format includes:
    - chrom, start, end: Genomic coordinates
    - name: CNV type and copy number (e.g., "DEL_CN1", "DUP_CN3")
    - score: Posterior probability * 1000

    Args:
        result: AnalysisResult containing CNV segments
        output_path: Path to output BED file
        cnv_only: If True, only write CNV calls
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Writing BED to {output_path}")

    with open(output_path, 'w') as f:
        # Write header
        f.write("#chrom\tstart\tend\tname\tscore\tstrand\n")

        count = 0
        for segment in result.segments:
            if cnv_only and not segment.is_cnv:
                continue

            sv_type = _get_sv_type(segment)
            name = f"{sv_type}_CN{segment.absolute_cn}"
            score = int(segment.probability * 1000)

            record = BEDRecord(
                chrom=segment.chrom,
                start=segment.start,
                end=segment.end,
                name=name,
                score=score
            )
            f.write(record.to_line() + "\n")
            count += 1

    logger.info(f"Wrote {count} records to BED")


def write_segments_tsv(
    result: AnalysisResult,
    output_path: Union[str, Path]
) -> None:
    """Write per-window results to TSV format (HSLM output format).

    This format includes one row per window, matching the original
    EXCAVATOR2 output format where each exon/window has:
    - Log2R: The raw log2 ratio for that specific window
    - SegMean: The mean log2 ratio of the segment the window belongs to

    Columns:
    - Chromosome, Position, Start, End, Log2R, SegMean, Class,
      CN, AbsoluteCN, Probability

    Args:
        result: AnalysisResult containing window_results
        output_path: Path to output TSV file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Writing segments TSV to {output_path}")

    with open(output_path, 'w') as f:
        # Write header
        headers = [
            "Chromosome", "Position", "Start", "End", "Log2R", "SegMean",
            "Class", "CN", "AbsoluteCN", "Probability"
        ]
        f.write("\t".join(headers) + "\n")

        # Use window_results if available, otherwise fall back to segments
        if result.window_results:
            for wr in result.window_results:
                # Format segment_mean - use "NA" if NaN (segmentation failed)
                if np.isnan(wr.segment_mean):
                    seg_mean_str = "NA"
                else:
                    seg_mean_str = f"{wr.segment_mean:.6f}"

                values = [
                    wr.chrom,
                    str(wr.position),
                    str(wr.start),
                    str(wr.end),
                    f"{wr.log2_ratio:.6f}",
                    seg_mean_str,
                    wr.region_class,
                    str(wr.cn_call),
                    str(wr.absolute_cn),
                    f"{wr.probability:.6f}"
                ]
                f.write("\t".join(values) + "\n")

            logger.info(f"Wrote {len(result.window_results)} window records to TSV")
        else:
            # Fallback: write segment-level data (legacy behavior)
            logger.warning("No window_results available, falling back to segment-level output")
            for segment in result.segments:
                position = (segment.start + segment.end) // 2

                values = [
                    segment.chrom,
                    str(position),
                    str(segment.start),
                    str(segment.end),
                    f"{segment.segment_mean:.6f}",
                    f"{segment.segment_mean:.6f}",
                    segment.region_class,
                    str(segment.cn_call),
                    str(segment.absolute_cn),
                    f"{segment.probability:.6f}"
                ]
                f.write("\t".join(values) + "\n")

            logger.info(f"Wrote {len(result.segments)} segment records to TSV")


def write_fastcall_bed(
    result: AnalysisResult,
    output_path: Union[str, Path]
) -> None:
    """Write FastCall results in extended BED format.

    This format matches the original EXCAVATOR2 FastCall output:
    Chromosome, Start, End, Segment(Log2R), CNF, CN, Call, ProbCall

    Args:
        result: AnalysisResult containing segments
        output_path: Path to output file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Writing FastCall BED to {output_path}")

    with open(output_path, 'w') as f:
        # Write header
        headers = [
            "Chromosome", "Start", "End", "Segment", "CNF", "CN", "Call", "ProbCall"
        ]
        f.write("\t".join(headers) + "\n")

        for segment in result.segments:
            # CNF is calculated copy number (before rounding)
            # For simplicity, we use absolute_cn directly
            cnf = segment.absolute_cn

            values = [
                segment.chrom,
                str(segment.start),
                str(segment.end),
                f"{segment.segment_mean:.6f}",
                str(cnf),
                str(segment.absolute_cn),
                str(segment.cn_call),
                f"{segment.probability:.6f}"
            ]
            f.write("\t".join(values) + "\n")

    logger.info(f"Wrote {len(result.segments)} records to FastCall BED")
