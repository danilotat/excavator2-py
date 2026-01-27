"""
Unit tests for EXCAVATOR2 target initialization module.
"""

import pytest
import numpy as np
import tempfile
from pathlib import Path

# Import IO modules
from excavator2.io.bed import (
    TargetRegion,
    ChromosomeInfo,
    GapRegion,
    load_target_bed,
    load_chromosome_coordinates,
    load_gap_file,
    regions_to_bed,
    merge_overlapping_regions,
    get_chromosome_regions,
    normalize_chromosome_name,
)

# Import target modules
from excavator2.target.filter import (
    AnalysisWindow,
    FilteredTargetResult,
    create_analysis_windows,
    filter_gap_overlapping_windows,
    create_filtered_target,
    renumber_windows,
)
from excavator2.target.gc_content import (
    get_window_gc_stats,
)
from excavator2.target.mappability import (
    get_mappability_stats,
    filter_low_mappability_windows,
)
from excavator2.target.init import (
    TargetData,
    save_target_data,
    load_target_data,
)


# =============================================================================
# Tests for io/bed.py
# =============================================================================

class TestTargetRegion:
    """Tests for TargetRegion dataclass."""

    def test_basic_region(self):
        region = TargetRegion("chr1", 1000, 2000)
        assert region.chrom == "chr1"
        assert region.start == 1000
        assert region.end == 2000
        assert region.length == 1000
        assert region.midpoint == 1500

    def test_region_with_name(self):
        region = TargetRegion("chr1", 1000, 2000, name="gene1")
        assert region.name == "gene1"

    def test_region_with_all_fields(self):
        region = TargetRegion("chr1", 1000, 2000, "gene1", 100.0, "+")
        assert region.score == 100.0
        assert region.strand == "+"


class TestChromosomeInfo:
    """Tests for ChromosomeInfo dataclass."""

    def test_basic_chrom(self):
        chrom = ChromosomeInfo("chr1", 0, 248956422)
        assert chrom.name == "chr1"
        assert chrom.start == 0
        assert chrom.end == 248956422
        assert chrom.length == 248956422


class TestGapRegion:
    """Tests for GapRegion dataclass."""

    def test_basic_gap(self):
        gap = GapRegion("chr1", 121535434, 124535434, "centromere")
        assert gap.chrom == "chr1"
        assert gap.length == 3000000
        assert gap.name == "centromere"


class TestLoadTargetBed:
    """Tests for load_target_bed function."""

    def test_load_simple_bed(self, tmp_path):
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("chr1\t1000\t2000\nghr1\t3000\t4000\n")

        regions = load_target_bed(bed_file)
        assert len(regions) == 2
        assert regions[0].chrom == "chr1"
        assert regions[0].start == 1000
        assert regions[0].end == 2000

    def test_load_bed_with_names(self, tmp_path):
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("chr1\t1000\t2000\tgene1\nchr1\t3000\t4000\tgene2\n")

        regions = load_target_bed(bed_file)
        assert regions[0].name == "gene1"
        assert regions[1].name == "gene2"

    def test_skip_header(self, tmp_path):
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("# header\nchr1\t1000\t2000\n")

        regions = load_target_bed(bed_file)
        assert len(regions) == 1

    def test_skip_track_lines(self, tmp_path):
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("track name=test\nbrowser position chr1\nchr1\t1000\t2000\n")

        regions = load_target_bed(bed_file)
        assert len(regions) == 1

    def test_filter_by_chromosome(self, tmp_path):
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("chr1\t1000\t2000\nchr2\t3000\t4000\nchr3\t5000\t6000\n")

        regions = load_target_bed(bed_file, chromosomes=["chr1", "chr2"])
        assert len(regions) == 2

    def test_filter_min_length(self, tmp_path):
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("chr1\t1000\t1050\nchr1\t2000\t3000\n")

        regions = load_target_bed(bed_file, min_length=100)
        assert len(regions) == 1
        assert regions[0].start == 2000


class TestLoadChromosomeCoordinates:
    """Tests for load_chromosome_coordinates function."""

    def test_load_coordinates(self, tmp_path):
        coord_file = tmp_path / "chromosomes.txt"
        coord_file.write_text("CHR\tSTART\tEND\nchr1\t1\t248956422\nchr2\t1\t242193529\n")

        chroms = load_chromosome_coordinates(coord_file)
        assert len(chroms) == 2
        assert "chr1" in chroms
        assert chroms["chr1"].end == 248956422


class TestLoadGapFile:
    """Tests for load_gap_file function."""

    def test_load_gaps(self, tmp_path):
        gap_file = tmp_path / "gaps.txt"
        gap_file.write_text("bin\tchrom\tchromStart\tchromEnd\tix\tn\tsize\ttype\tbridge\n"
                           "1\tchr1\t0\t10000\t1\tN\t10000\ttelomere\tno\n"
                           "2\tchr1\t121535434\t124535434\t2\tN\t3000000\tcentromere\tno\n")

        gaps = load_gap_file(gap_file)
        assert len(gaps) == 2
        assert gaps[0].name == "telomere"
        assert gaps[1].name == "centromere"

    def test_exclude_alt_chromosomes(self, tmp_path):
        gap_file = tmp_path / "gaps.txt"
        gap_file.write_text("bin\tchrom\tchromStart\tchromEnd\tix\tn\tsize\ttype\tbridge\n"
                           "1\tchr1\t0\t10000\t1\tN\t10000\ttelomere\tno\n"
                           "2\tchr1_random\t0\t5000\t2\tN\t5000\tother\tno\n")

        gaps = load_gap_file(gap_file, exclude_alt=True)
        assert len(gaps) == 1


class TestRegionsToBed:
    """Tests for regions_to_bed function."""

    def test_write_bed(self, tmp_path):
        regions = [
            TargetRegion("chr1", 1000, 2000, "gene1"),
            TargetRegion("chr1", 3000, 4000, "gene2"),
        ]
        output_file = tmp_path / "output.bed"

        regions_to_bed(regions, output_file)

        content = output_file.read_text()
        lines = content.strip().split("\n")
        assert len(lines) == 2
        assert lines[0].startswith("chr1\t1000\t2000")


class TestMergeOverlappingRegions:
    """Tests for merge_overlapping_regions function."""

    def test_merge_overlapping(self):
        regions = [
            TargetRegion("chr1", 1000, 2000),
            TargetRegion("chr1", 1500, 2500),
            TargetRegion("chr1", 4000, 5000),
        ]

        merged = merge_overlapping_regions(regions)
        assert len(merged) == 2
        assert merged[0].end == 2500

    def test_merge_adjacent(self):
        regions = [
            TargetRegion("chr1", 1000, 2000),
            TargetRegion("chr1", 2000, 3000),
        ]

        merged = merge_overlapping_regions(regions)
        assert len(merged) == 1
        assert merged[0].end == 3000

    def test_no_merge_different_chroms(self):
        regions = [
            TargetRegion("chr1", 1000, 2000),
            TargetRegion("chr2", 1000, 2000),
        ]

        merged = merge_overlapping_regions(regions)
        assert len(merged) == 2


class TestNormalizeChromosomeName:
    """Tests for normalize_chromosome_name function."""

    def test_add_chr_prefix(self):
        assert normalize_chromosome_name("1", use_chr_prefix=True) == "chr1"
        assert normalize_chromosome_name("chr1", use_chr_prefix=True) == "chr1"

    def test_remove_chr_prefix(self):
        assert normalize_chromosome_name("chr1", use_chr_prefix=False) == "1"
        assert normalize_chromosome_name("1", use_chr_prefix=False) == "1"


# =============================================================================
# Tests for target/filter.py
# =============================================================================

class TestAnalysisWindow:
    """Tests for AnalysisWindow dataclass."""

    def test_basic_window(self):
        window = AnalysisWindow("chr1", 1000, 2000, "a1", "IN")
        assert window.chrom == "chr1"
        assert window.start == 1000
        assert window.end == 2000
        assert window.window_id == "a1"
        assert window.region_class == "IN"
        assert window.length == 1000
        assert window.position == 1500

    def test_custom_position(self):
        window = AnalysisWindow("chr1", 1000, 2000, "a1", "OUT", position=1200)
        assert window.position == 1200


class TestFilteredTargetResult:
    """Tests for FilteredTargetResult dataclass."""

    def test_basic_result(self):
        windows = [
            AnalysisWindow("chr1", 1000, 2000, "a1", "IN"),
            AnalysisWindow("chr1", 3000, 4000, "a2", "OUT"),
        ]
        result = FilteredTargetResult(
            windows=windows,
            chromosomes=["chr1"],
            n_in_target=1,
            n_out_target=1,
            n_filtered=0,
            target_name="test"
        )
        assert result.n_windows == 2

    def test_get_chromosome_windows(self):
        windows = [
            AnalysisWindow("chr1", 1000, 2000, "a1", "IN"),
            AnalysisWindow("chr2", 3000, 4000, "a2", "IN"),
        ]
        result = FilteredTargetResult(
            windows=windows,
            chromosomes=["chr1", "chr2"],
            n_in_target=2,
            n_out_target=0,
            n_filtered=0,
            target_name="test"
        )

        chr1_windows = result.get_chromosome_windows("chr1")
        assert len(chr1_windows) == 1
        assert chr1_windows[0].chrom == "chr1"


class TestCreateAnalysisWindows:
    """Tests for create_analysis_windows function."""

    def test_in_target_windows(self):
        regions = [
            TargetRegion("chr1", 1000, 2000),
            TargetRegion("chr1", 3000, 4000),
        ]
        chrom_info = {"chr1": ChromosomeInfo("chr1", 0, 10000)}

        windows = create_analysis_windows(regions, chrom_info, window_size=1000)

        # Should have at least the 2 IN-target windows
        in_windows = [w for w in windows if w.region_class == "IN"]
        assert len(in_windows) == 2

    def test_out_target_windows(self):
        regions = [
            TargetRegion("chr1", 1000, 2000),
            TargetRegion("chr1", 50000, 51000),  # Large gap between regions
        ]
        chrom_info = {"chr1": ChromosomeInfo("chr1", 0, 100000)}

        windows = create_analysis_windows(regions, chrom_info, window_size=1000, flank=200)

        # Should have OUT-target windows in the gap
        out_windows = [w for w in windows if w.region_class == "OUT"]
        assert len(out_windows) > 0


class TestFilterGapOverlappingWindows:
    """Tests for filter_gap_overlapping_windows function."""

    def test_filter_overlapping(self):
        windows = [
            AnalysisWindow("chr1", 1000, 2000, "a1", "IN"),
            AnalysisWindow("chr1", 5000, 6000, "a2", "IN"),  # Overlaps gap
            AnalysisWindow("chr1", 10000, 11000, "a3", "IN"),
        ]
        gaps = [GapRegion("chr1", 4000, 7000, "gap")]

        filtered, n_filtered = filter_gap_overlapping_windows(windows, gaps)
        assert len(filtered) == 2
        assert n_filtered == 1

    def test_no_gaps(self):
        windows = [
            AnalysisWindow("chr1", 1000, 2000, "a1", "IN"),
        ]
        gaps = []

        filtered, n_filtered = filter_gap_overlapping_windows(windows, gaps)
        assert len(filtered) == 1
        assert n_filtered == 0


class TestRenumberWindows:
    """Tests for renumber_windows function."""

    def test_renumber(self):
        windows = [
            AnalysisWindow("chr1", 1000, 2000, "x1", "IN"),
            AnalysisWindow("chr1", 3000, 4000, "x5", "OUT"),
        ]

        renumbered = renumber_windows(windows)
        assert renumbered[0].window_id == "a1"
        assert renumbered[1].window_id == "a2"


# =============================================================================
# Tests for target/gc_content.py
# =============================================================================

class TestGetWindowGCStats:
    """Tests for get_window_gc_stats function."""

    def test_basic_stats(self):
        gc_content = np.array([0.4, 0.5, 0.6, 0.45, 0.55])
        windows = [
            AnalysisWindow("chr1", i*1000, (i+1)*1000, f"a{i}", "IN")
            for i in range(5)
        ]

        stats = get_window_gc_stats(gc_content, windows)
        assert 'mean' in stats
        assert 'median' in stats
        assert 'std' in stats
        assert 0.4 <= stats['mean'] <= 0.6

    def test_empty_gc(self):
        gc_content = np.array([])
        windows = []

        stats = get_window_gc_stats(gc_content, windows)
        assert stats['mean'] == 0.0

    def test_in_out_target_stats(self):
        gc_content = np.array([0.3, 0.4, 0.5, 0.6])
        windows = [
            AnalysisWindow("chr1", 0, 1000, "a1", "IN"),
            AnalysisWindow("chr1", 1000, 2000, "a2", "IN"),
            AnalysisWindow("chr1", 2000, 3000, "a3", "OUT"),
            AnalysisWindow("chr1", 3000, 4000, "a4", "OUT"),
        ]

        stats = get_window_gc_stats(gc_content, windows)
        assert stats['in_target_mean'] == 0.35
        assert stats['out_target_mean'] == 0.55


# =============================================================================
# Tests for target/mappability.py
# =============================================================================

class TestGetMappabilityStats:
    """Tests for get_mappability_stats function."""

    def test_basic_stats(self):
        mappability = np.array([0.8, 0.9, 1.0, 0.95, 0.85])
        windows = [
            AnalysisWindow("chr1", i*1000, (i+1)*1000, f"a{i}", "IN")
            for i in range(5)
        ]

        stats = get_mappability_stats(mappability, windows)
        assert 'mean' in stats
        assert 'low_mappability_count' in stats

    def test_low_mappability_count(self):
        mappability = np.array([0.3, 0.4, 0.9, 1.0])
        windows = [
            AnalysisWindow("chr1", i*1000, (i+1)*1000, f"a{i}", "IN")
            for i in range(4)
        ]

        stats = get_mappability_stats(mappability, windows)
        assert stats['low_mappability_count'] == 2  # < 0.5 threshold


class TestFilterLowMappabilityWindows:
    """Tests for filter_low_mappability_windows function."""

    def test_filter_low_map(self):
        windows = [
            AnalysisWindow("chr1", 0, 1000, "a1", "IN"),
            AnalysisWindow("chr1", 1000, 2000, "a2", "IN"),
            AnalysisWindow("chr1", 2000, 3000, "a3", "IN"),
        ]
        mappability = np.array([0.05, 0.8, 0.9])

        filtered_windows, filtered_map, n_filtered = filter_low_mappability_windows(
            windows, mappability, min_mappability=0.1
        )
        assert len(filtered_windows) == 2
        assert n_filtered == 1


# =============================================================================
# Tests for target/init.py
# =============================================================================

class TestTargetData:
    """Tests for TargetData dataclass."""

    def test_basic_target_data(self):
        windows = [
            AnalysisWindow("chr1", 0, 1000, "a1", "IN"),
            AnalysisWindow("chr1", 2000, 3000, "a2", "OUT"),
            AnalysisWindow("chr2", 0, 1000, "a3", "IN"),
        ]
        gc_content = np.array([0.4, 0.5, 0.45])
        mappability = np.array([0.9, 0.8, 0.95])

        target = TargetData(
            windows=windows,
            gc_content=gc_content,
            mappability=mappability,
            chromosomes=["chr1", "chr2"],
            target_name="test",
            assembly="hg38",
            window_size=1000
        )

        assert target.n_windows == 3
        assert target.n_in_target == 2
        assert target.n_out_target == 1

    def test_get_chromosome_data(self):
        windows = [
            AnalysisWindow("chr1", 0, 1000, "a1", "IN"),
            AnalysisWindow("chr2", 0, 1000, "a2", "IN"),
        ]
        gc_content = np.array([0.4, 0.5])
        mappability = np.array([0.9, 0.8])

        target = TargetData(
            windows=windows,
            gc_content=gc_content,
            mappability=mappability,
            chromosomes=["chr1", "chr2"],
            target_name="test",
            assembly="hg38",
            window_size=1000
        )

        chr1_data = target.get_chromosome_data("chr1")
        assert chr1_data['n_windows'] == 1
        assert chr1_data['gc_content'][0] == 0.4


class TestTargetDataHDF5:
    """Tests for save/load target data HDF5 functions."""

    def test_save_load_roundtrip(self, tmp_path):
        windows = [
            AnalysisWindow("chr1", 0, 1000, "a1", "IN"),
            AnalysisWindow("chr1", 2000, 3000, "a2", "OUT"),
        ]
        gc_content = np.array([0.4, 0.5])
        mappability = np.array([0.9, 0.8])

        original = TargetData(
            windows=windows,
            gc_content=gc_content,
            mappability=mappability,
            chromosomes=["chr1"],
            target_name="test",
            assembly="hg38",
            window_size=1000,
            metadata={'creation_date': '2024-01-01'}
        )

        # Save
        output_file = tmp_path / "target.h5"
        save_target_data(original, output_file)
        assert output_file.exists()

        # Load
        loaded = load_target_data(output_file)

        assert loaded.target_name == original.target_name
        assert loaded.assembly == original.assembly
        assert loaded.window_size == original.window_size
        assert loaded.n_windows == original.n_windows
        assert len(loaded.chromosomes) == len(original.chromosomes)
        np.testing.assert_array_almost_equal(loaded.gc_content, original.gc_content)
        np.testing.assert_array_almost_equal(loaded.mappability, original.mappability)

    def test_save_empty_target(self, tmp_path):
        target = TargetData(
            windows=[],
            gc_content=np.array([]),
            mappability=np.array([]),
            chromosomes=[],
            target_name="empty",
            assembly="hg38",
            window_size=1000
        )

        output_file = tmp_path / "empty.h5"
        save_target_data(target, output_file)
        assert output_file.exists()

        loaded = load_target_data(output_file)
        assert loaded.n_windows == 0
