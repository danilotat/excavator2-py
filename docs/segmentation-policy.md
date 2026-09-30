# Segmentation policy

This replaces the legacy wrapper policy for H04 (#12), A22 (#10), A24 (#16),
A23 (#15), A04 (#19), and A05 (#3). There is no legacy mode. The frozen numerical
kernel fixtures still verify transitions and Viterbi independently of wrapper
policy. The original R/Fortran sources and fixture files remain historical evidence.

## Noise and state space

Analysis estimates parameters independently within each chromosome arm. Emission
noise is `median(abs(diff(log2_ratio))) / (sqrt(2) * Phi^-1(0.75))`, with a
0.001 log2-ratio floor. With IN/OUT metadata, each class uses differences between
successive observations of that class in the arm. Classes with fewer than three
observations use the pooled arm estimate. Constant profiles and single-window
arms have finite noise and retain their observed medians.

This estimator targets independent Gaussian measurement errors with sparse level
changes. Unlike variance of profile levels, distant chromosome CNV burden cannot
change the local arm's emission noise, prior, or state space. It does not model
per-window coverage, autocorrelation, or frequent alternating biological shifts;
class-specific precision is not a claim of calibrated posterior uncertainty.

`Omega` now specifies the Gaussian reset-prior variance in squared log2-ratio
units (`smu = sqrt(Omega)`, center zero); its accepted range remains `(0, 1)` and
default remains 0.1. It no longer divides global profile variance into signal and
noise. This explicit prior decouples regularization from CNV burden and avoids
collapsing the reset distribution when measured noise approaches zero. Existing
parameter files with tuned Omega values need reassessment under this interpretation.

The candidate states are zero plus occupied bins of the arm's observed ratios,
with spacing `max(0.001, min(0.05, minimum_class_noise / 2))`. Every observation is
represented within half a bin; the amplitude range is not clipped. Unoccupied
bins are omitted to reduce kernel cost. The discrete Gaussian is normalized on
these candidate states. State selection is data dependent; it is not a continuous
state model. Final segment values remain medians of original observations.
Kernel runtime is O(windows × states²), with O(windows × states) traceback storage;
full diagnostic traces additionally store the transition matrices.

## Arm prior and filtering

The initial state uses the same normalized Gaussian distribution as transition
resets. The resulting transition model is reversible: reversing observations,
class precision, and inter-window distances preserves path probabilities. An
exactly tied optimum can still be resolved differently by Viterbi's deterministic
first-state tie rule; tests assert reversal for unique optima, not arbitrary ties.

`minExons` means **at least that many IN windows**, not all windows and not distinct
exon IDs (which are absent from prepared artifacts). Exactly-threshold segments
survive. The low-level `segment` API counts all observations unless a boolean
support vector is supplied; analysis always supplies the IN mask.

Among boundaries touching an insufficient-support segment, remove those with
the smallest absolute difference between neighboring segment medians, then
recompute support and medians. Equally good boundaries are removed together;
this avoids choosing a coordinate direction even for symmetric ties. A tie can
merge three or more neighboring segments. This is a median-distance merge policy,
not a significance test. Filtering never crosses arm boundaries and never removes
the first or final endpoint. An entire arm below the support threshold is retained
because it has no legal neighbor; its median is never replaced with zero.

## Segment identities and outputs

Every filtered segment receives a unique ID across arms and chromosomes. Analysis
passes these IDs to summarization and stores them in `checkpoints.npz`; table and
VCF generation consume those explicit intervals. Identical medians therefore do
not erase inferred boundaries, centromeres, or chromosome boundaries. The public
`segment_profile` return tuple now includes IDs as its fourth member, and
`summarize` requires IDs. Output table columns are unchanged.

The analysis manifest records the segmentation policy and support units. Existing
rejections for malformed centromere geometry remain in place.

## Validation and changed calls

`tests/test_segmentation_policy.py` covers positive/negative 1.2→2.4 plateaus,
0.31→0.34 steps, unchanged local-event detection under distant ±1/±2 burden,
weak/strong terminal reversal with unequal distances, neutral-noise false positives,
IN/OUT precision, native/reference agreement, first/middle/last/sole filtering,
exact support thresholds, and equal medians across arms/chromosomes.

The frozen non-normal profiles were also evaluated without changing their input
ratios. Both old and new segmentations were classified by the same current
FastCall implementation, isolating segmentation's effect from the other recent
normalization and classifier corrections:

| Profile | Segments, old→new | Called segments, old→new | Changed window labels |
|---|---:|---:|---:|
| paired/Test1 | 4→12 | 4→8 | 80/240 |
| paired/Test2 | 4→12 | 4→8 | 120/240 |
| pooling/Test1 | 4→12 | 4→8 | 80/240 |
| pooling/Test2 | 4→12 | 4→8 | 120/240 |

Full levels, boundaries, labels, and changed row indices are in
[segmentation-report.json](segmentation-report.json). Reproduce with:

```sh
python tools/validation/segmentation_report.py > docs/segmentation-report.json
pytest tests/test_segmentation_policy.py tests/test_hslm.py tests/test_analysis.py
```

These are synthetic/frozen regression profiles, not a held-out clinical cohort.
The paired and pooled profiles are not independent evidence. No claim is made
that every changed call is biologically correct or that sensitivity, specificity,
or probability calibration is established for real non-normal samples. Such cohort
qualification remains outstanding; the report makes the intentional changes visible.
