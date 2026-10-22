"""Merging of adjacent same-label segments into reported events (#76)."""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from excavator2.fastcall import call_probabilities, fastcall_events, segment_standard_errors

NOISE = {"IN": 0.1, "OUT": 0.2}


def events(levels, sizes, labels, arms=None, k=2.0, classes=None, cellularity=1.0):
    values = np.repeat(np.asarray(levels, dtype=float), sizes)
    ends = np.cumsum(sizes)
    starts = ends - sizes
    classes = np.array(["IN"] * len(values)) if classes is None else classes
    arms = np.array(["chr1:p"] * len(levels)) if arms is None else np.asarray(arms)
    return fastcall_events(
        values,
        classes,
        starts,
        ends,
        np.asarray(levels, dtype=float),
        np.asarray(labels),
        arms,
        NOISE,
        k=k,
        cellularity=cellularity,
    )


def test_consistent_gains_merge_and_are_rescored():
    rng = np.random.default_rng(3)
    # Two weak gains near the threshold, with alternating IN/OUT noise.
    noise = np.tile([0.1, 0.2], 30) * rng.standard_normal(60)
    values = np.repeat([0.40, 0.42], 30) + noise
    classes = np.array(["IN", "OUT"] * 30)
    starts, ends = np.array([0, 30]), np.array([30, 60])
    segments = np.array([np.median(values[:30]), np.median(values[30:])])
    result = fastcall_events(values, classes, starts, ends, segments, [1, 1], ["chr1:p"] * 2, NOISE)
    assert_array_equal(result.ids, [0, 0])
    assert_array_equal([result.starts, result.ends, result.labels], [[0], [60], [1]])
    assert_allclose(result.mixture, np.median(values))
    error = segment_standard_errors(classes, [0], [60], NOISE)
    assert_allclose(result.probabilities, call_probabilities([np.median(values)], error))
    pieces = call_probabilities(segments, segment_standard_errors(classes, starts, ends, NOISE))
    # More windows give a smaller standard error: the weaker piece gains certainty.
    assert pieces[:, 3].min() < 0.7 < 0.95 < result.probabilities[0, 3]


def test_significant_steps_labels_and_arms_stay_separate():
    # A within-class step, supported by many windows, is a real level change.
    assert_array_equal(events([0.45, 0.62], [200, 200], [1, 1]).ids, [0, 1])
    # Different calls never merge, however close their medians.
    assert_array_equal(events([0.34, 0.36], [5, 5], [0, 1]).ids, [0, 1])
    # Identical medians do not merge across a centromere.
    assert_array_equal(events([0.5, 0.5], [5, 5], [1, 1], ["chr1:p", "chr1:q"]).ids, [0, 1])
    # MergeK 0 merges only identical medians.
    assert_array_equal(events([0.5, 0.5, 0.51], [5, 5, 5], [1, 1, 1], k=0).ids, [0, 0, 1])


@pytest.mark.parametrize("reverse", [False, True])
def test_most_consistent_pair_merges_first_and_statistics_are_recomputed(reverse):
    # B joins A first (z = 0). B alone agrees with C (z = 1.4), but A+B does not (z = 4.6).
    levels, sizes = np.array([0.5, 0.5, 0.56]), np.array([1000, 10, 100])
    expected = [0, 0, 1]
    if reverse:
        levels, sizes, expected = levels[::-1], sizes[::-1], [0, 1, 1]
    result = events(levels, sizes, [1, 1, 1])
    assert_array_equal(result.ids, expected)


def test_equal_ties_merge_together_and_unmerged_values_are_kept():
    result = events([0.5, 0.52, 0.5, 0.9 - 1e-9, -0.6], [20, 20, 20, 3, 5], [1, 1, 1, 1, -1])
    assert_array_equal(result.ids, [0, 0, 0, 1, 2])
    assert_allclose(result.mixture[1:], [0.9 - 1e-9, -0.6], rtol=0, atol=0)


def test_cellularity_applies_to_event_levels():
    result = events([0.3, 0.3], [10, 10], [1, 1], cellularity=0.5)
    assert_allclose(result.corrected, np.log2(2 * 2**0.3 - 1))
    error = np.sqrt(np.pi / 2) * 0.1 / np.sqrt(20)
    assert_allclose(result.probabilities, call_probabilities([0.3], [error], cellularity=0.5))


def test_invalid_merge_threshold_is_rejected():
    with pytest.raises(ValueError, match="MergeK"):
        events([0.5, 0.5], [5, 5], [1, 1], k=-1)
