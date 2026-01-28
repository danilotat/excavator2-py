"""
Shared utilities for EXCAVATOR2 plotting.

Contains color schemes, styling constants, and helper functions
used across all plotting modules.
"""

from typing import Dict, Tuple, List
import numpy as np

# CNV Call colors (matching original EXCAVATOR2)
# Using colorblind-friendly palette
CNV_COLORS: Dict[int, Tuple[float, float, float]] = {
    -2: (213/255, 94/255, 0/255),     # 2-copy deletion (dark orange)
    -1: (230/255, 159/255, 0/255),    # 1-copy deletion (light orange)
    0: (0.7, 0.7, 0.7),                # Normal (grey)
    1: (86/255, 180/255, 233/255),    # 1-copy gain (light blue)
    2: (0/255, 114/255, 178/255),     # N-copy amplification (dark blue)
}

# Labels for CNV states
CNV_LABELS: Dict[int, str] = {
    -2: "2-copy DEL",
    -1: "1-copy DEL",
    0: "Normal",
    1: "1-copy AMP",
    2: "N-copy AMP",
}

# Short labels for legends
CNV_SHORT_LABELS: Dict[int, str] = {
    -2: "2-DEL",
    -1: "DEL",
    0: "Normal",
    1: "AMP",
    2: "2-AMP",
}

# Region class colors
REGION_COLORS: Dict[str, str] = {
    "IN": "#1f77b4",       # Blue for in-target
    "OUT": "#87ceeb",      # Light blue for off-target
}

# Plot styling constants
FIGURE_DPI = 150
DEFAULT_FIGSIZE = (15, 10)
SCATTER_SIZE = 3
SCATTER_ALPHA = 0.6
LINE_WIDTH = 1.5
SEGMENT_LINE_WIDTH = 2.5
SEGMENT_COLOR = "#d62728"  # Red for segment line

# Axis limits
Y_LIM_LOG2 = (-3.0, 3.0)

# Chromosome order for genome-wide plots (GRCh38)
CHROMOSOME_ORDER = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
]

# Chromosome lengths (GRCh38/hg38) for genome-wide plots
CHROMOSOME_LENGTHS_HG38: Dict[str, int] = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}


def get_cnv_color(cn_call: int) -> Tuple[float, float, float]:
    """Get color for a CNV call.

    Args:
        cn_call: CNV call value (-2, -1, 0, 1, 2)

    Returns:
        RGB color tuple (0-1 range)
    """
    # Clamp to valid range
    cn_call = max(-2, min(2, cn_call))
    return CNV_COLORS.get(cn_call, CNV_COLORS[0])


def get_cnv_label(cn_call: int, short: bool = False) -> str:
    """Get label for a CNV call.

    Args:
        cn_call: CNV call value
        short: Use short labels (for legends)

    Returns:
        String label
    """
    cn_call = max(-2, min(2, cn_call))
    labels = CNV_SHORT_LABELS if short else CNV_LABELS
    return labels.get(cn_call, "Unknown")


def compute_bin_statistics(
    values: np.ndarray,
    counts: np.ndarray,
    bin_size: float,
    min_count: int = 10
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute binned statistics for bias plots.

    Calculates median and quartiles for counts in each value bin.

    Args:
        values: Feature values to bin by (e.g., GC content)
        counts: Read counts
        bin_size: Size of each bin
        min_count: Minimum windows in bin to include

    Returns:
        Tuple of (bin_centers, medians, q1, q3)
    """
    if len(values) == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])

    # Create bins
    val_min = np.floor(np.min(values) / bin_size) * bin_size
    val_max = np.ceil(np.max(values) / bin_size) * bin_size
    bins = np.arange(val_min, val_max + bin_size, bin_size)

    bin_centers = []
    medians = []
    q1_values = []
    q3_values = []

    for i in range(len(bins) - 1):
        mask = (values >= bins[i]) & (values < bins[i + 1])
        bin_counts = counts[mask]

        if len(bin_counts) >= min_count:
            # Filter out zeros for statistics
            nonzero = bin_counts[bin_counts > 0]
            if len(nonzero) >= min_count // 2:
                bin_centers.append(bins[i] + bin_size / 2)
                medians.append(np.median(nonzero))
                q1_values.append(np.percentile(nonzero, 25))
                q3_values.append(np.percentile(nonzero, 75))

    return (
        np.array(bin_centers),
        np.array(medians),
        np.array(q1_values),
        np.array(q3_values)
    )


def normalize_chromosome_name(chrom: str) -> str:
    """Normalize chromosome name to chr-prefixed format.

    Args:
        chrom: Chromosome name (e.g., "1", "chr1", "Chr1")

    Returns:
        Normalized name (e.g., "chr1")
    """
    chrom = str(chrom).lower()
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    return chrom


def sort_chromosomes(chromosomes: List[str]) -> List[str]:
    """Sort chromosome names in standard order.

    Args:
        chromosomes: List of chromosome names

    Returns:
        Sorted list
    """
    def sort_key(chrom):
        norm = normalize_chromosome_name(chrom)
        if norm in CHROMOSOME_ORDER:
            return (0, CHROMOSOME_ORDER.index(norm))
        # Put unknowns at the end
        return (1, chrom)

    return sorted(chromosomes, key=sort_key)


def format_position(pos: int) -> str:
    """Format genomic position for display.

    Args:
        pos: Position in base pairs

    Returns:
        Formatted string (e.g., "12.5 Mb")
    """
    if pos >= 1_000_000:
        return f"{pos / 1_000_000:.1f} Mb"
    elif pos >= 1_000:
        return f"{pos / 1_000:.1f} kb"
    else:
        return f"{pos} bp"
