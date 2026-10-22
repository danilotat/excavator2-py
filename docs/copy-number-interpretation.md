# Copy-number interpretation and remaining model limits

This resolves the reporting contract for A18 ([#6](https://github.com/danilotat/excavator2-py/issues/6))
and A25 ([#11](https://github.com/danilotat/excavator2-py/issues/11)), assesses
P15 ([#18](https://github.com/danilotat/excavator2-py/issues/18)), and addresses
H07 ([#17](https://github.com/danilotat/excavator2-py/issues/17)) in the current
Python pipeline. Historical R behavior is unchanged.

## A18: one corrected quantity for calling and reporting

For mixture segment log2 ratio m and cellularity p < 1, the corrected ratio is
`c = log2(max((2**m - (1-p))/p, 2**-5))`. At p = 1, c = m, with no floor.
FastCall receives c. The table exposes `MixtureSegment` and `CorrectedSegment`;
`DECNF = 2 * 2**c` and `DECN = rint(DECNF)` describe the corrected component.
Both VCF writers use exactly these values, without rounding the fraction before
rounding the integer. Thus at p = 0.5 the mixture fractions 1.5 and 2.5 produce
corrected values 1 and 3. HSLM rows and window `L2R` remain observed mixture
quantities. Checkpoints retain both segment arrays.

The floor for p < 1 produces DECNF = 0.0625 even for an impossible negative
inverted ratio. This is a numerical floor, not evidence for a small positive
biological copy count. Small purity amplifies ratio errors; neither correction
nor integer rounding supplies a confidence interval. `rint` uses ties to even.

## A25: explicit diploid-equivalent reporting, not absolute ploidy inference

The ambiguous `CNF`/`CN` fields are replaced by `DECNF`/`DECN` in tables and VCFs.
The manifest declares `absolute_copy_number: false`. VCF `GT` is missing (`.`),
since FastCall does not infer allele genotypes. There is no compatibility mode.
Consumers must update their field selection; historical output equality is no
longer an acceptance criterion.

Reference ploidy, sex and PAR boundaries are not inferred or validated by the
pipeline. DECNF/DECN are diploid-equivalent relative values, not absolute copy
numbers. For example, a twofold ratio relative to a haploid X reference produces
DECNF 4, not absolute CN 2. No user declaration of a diploid reference is required.
The contaminating component assumed by cellularity correction also affects how
its corrected values should be interpreted.

This addresses misleading absolute naming, not general ploidy-aware inference.
Adding absolute CN later requires region-specific reference and contaminant
ploidy, sample context, validated PAR handling, and an appropriate calling model.
The fixed FastCall means and class thresholds retain their diploid-relative meaning.

## Ploidy re-centering ([#75](https://github.com/danilotat/excavator2-py/issues/75))

Depth normalization and IN/OUT median centering set log2 ratio 0 to the copy
number of the median window, in paired and pooled designs alike; a matched normal
removes locus bias but not the genome-wide scale. FastCall reads that level as
diploid. In a tumor of ploidy 2.85 most of the genome is then shifted: CN 2 sits
near `log2(2/2.85) = -0.51` and is called a loss, while CN 3 sits near 0.

The optional FastCall `Ploidy` parameter (default 2, no shift) supplies the known
tumor ploidy P. After segmentation, each segment median m becomes
`m + log2((p*P + 2*(1-p)) / 2)`, where p is `Cellularity` and the contaminant is
diploid. Cellularity correction, FastCall labels, probabilities, `MixtureSegment`,
`CorrectedSegment` and DECNF use the shifted value. Segmentation does not, because
its reset prior is centered at zero; segment values are medians of observations,
so the shift is an exact translation. HSLM rows and window `L2R` stay
median-centered observations.

This is the manual re-centering used for CNVkit in published benchmarks, not
ploidy inference. It assumes the median window's copy number approximates the
mean ploidy; for a genome mostly at CN 3 with ploidy 2.85, the shift is off by
`log2(3/2.85) = 0.07`. DECNF remains diploid-equivalent and is absolute only to
the extent that P and this assumption hold.

## P15 assessment: control depth still affects arithmetic pooling

Analysis averages normalized counts across selected controls without a pseudocount,
and median-centers the resulting log ratios separately for IN and OUT. This is
the selected original automatic normalization method; no exposure file is used.

A deeper control has greater influence on local deviations in the pooled
reference. With one flat control and another carrying a twofold gain in ten of
100 windows, increasing only the second control's depth changes the local ratio
even after median centering. The regression test records this limitation rather
than claiming depth-invariant pooling. Control CNVs, quality weighting and outlier
handling remain unresolved modeling questions; no replacement pooling rule is
introduced. See [automatic normalization](analysis-normalization.md).

## H07: probabilities from segment measurement uncertainty

The FastCall fit sees one median per segment. Since A19, its truncated class
supports have disjoint interiors, so its posterior is 1 for the class containing
the median, irrespective of length or noise. Labels still come from that fit,
which equals the interval containing the median, and it is kept as a diagnostic
(`posterior`, `fastcall_trace`). Reported probabilities no longer use it.

`ProbCall`/`FCP` is the probability that the true segment level lies in the
called class interval. `ProbAltered`/`FPA` is the probability that it lies
outside the normal interval. Both come from each segment's own windows:

1. Per sample, estimate robust window noise separately for IN and OUT from
   adjacent same-class differences within segments:
   `sigma = median(|difference|) / 0.9539`. Steps between segments are excluded.
2. A segment median of `n_IN` and `n_OUT` windows has standard error
   `SE = sqrt(pi/2) * sqrt(n_IN * sigma_IN^2 + n_OUT * sigma_OUT^2) / n`.
3. With the observed mixture median normal around the true level and a flat
   prior, each class probability is a difference of normal CDFs at its edges.
   Corrected edges `-1.3, -d, u, 0.9` map to the mixture scale through the
   monotone cellularity correction. The outer classes are open-ended, so the
   five probabilities sum to one, and mixtures below the correction floor are
   deep losses.

Short or noisy segments, and segments near a threshold, therefore get lower
probabilities; fragmenting other parts of a profile does not change them. A
median just below 0.9 is a gain with `ProbCall` about 0.5 but `ProbAltered` 1.
`ProbCall` can be below another class's probability when the median lies near
the edge of a narrow interval; the label remains the interval containing it.
The manifest records `fastcall_probability_semantics:
segment-level-interval-probability` and per-sample `measurement_noise`;
checkpoints save the noise, standard errors and all five probabilities.

Windows are treated as independent. On two real paired WES runs of one cell line
(about 354,000 windows each), the lag-1 correlation of within-segment residuals
was about 0.05 (rank correlation about 0.03), at most about 0.09 for IN windows
under 1 kb apart. That reduces the effective window count by about 10%, so no
correction is applied. The noise level differed between runs (OUT 0.31 versus
0.57), which is why it is estimated per sample and class.

These probabilities describe random measurement error given the segmentation.
They exclude boundary uncertainty, purity uncertainty and systematic bias such
as a control CNV in the pool (P15), so they are not calibrated event confidence
and must not be used as empirical error rates. Calibration remains open: compare
calls between replicate runs, or against truth-labelled samples, by length,
depth, IN/OUT composition, noise and purity.

Regression evidence is in `tests/test_reporting_semantics.py`,
`tests/test_analysis.py` and `tests/test_signal_preservation.py`.

## Event merging ([#76](https://github.com/danilotat/excavator2-py/issues/76))

HSLM often splits one alteration into adjacent segments with the same FastCall
label, for example across a large gap between targets, where its distance prior
makes a boundary cheap. Analysis groups such segments into events, after the
merge pass of DRAGEN's
[CNV smoothing](https://help.dragen.illumina.com/dragen-v4.4/product-guide/dragen-v4.4/dragen-dna-pipeline/cnv-calling/additional-documentation/cnv-ascn#wgs-cnv-smoothing).

Two adjacent events in one chromosome arm, with the same label, merge when
`|m1 - m2| / sqrt(SE1^2 + SE2^2) <= MergeK`, using mixture medians and the H07
standard errors. The most consistent boundaries are removed first (equal ties
together), then the merged median and standard error are recomputed from all of
the event's windows, and the test repeats until no boundary qualifies. A merged
median lies between its parts' medians, so the label never changes. Events are
rescored like segments: DECNF and probabilities come from the merged median and
its standard error. An unmerged event keeps its segment's values exactly.

The test uses evidence, not distance. A bp merge distance would reflect target
layout on exome and panel data, and would block the gap-induced splits merging is
meant to undo. A within-class step supported by many windows, such as a subclonal
gain next to a clonal one, exceeds the threshold and stays two events. `MergeK`
is an optional FastCall parameter (default 2); 0 merges only identical medians.
DRAGEN's removal of negligible segments corresponds to `minExons`.

Merging is reporting only. HSLM rows, the segment-level `FastCallResults` rows and
the Window VCF keep every segment; `FastCallResults` gains an `EventId` column,
and the Region VCF has one record per called event. Checkpoints save event IDs,
window ranges, medians, labels and probabilities, and the manifest records
`event_policy`. Evidence is in `tests/test_events.py` and `tests/test_analysis.py`.
