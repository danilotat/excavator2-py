# Copy-number interpretation and remaining model limits

This resolves the reporting contract for A18 ([#6](https://github.com/danilotat/excavator2-py/issues/6))
and A25 ([#11](https://github.com/danilotat/excavator2-py/issues/11)), and assesses
P15 ([#18](https://github.com/danilotat/excavator2-py/issues/18)) and
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

## H07 assessment: probabilities are not calibrated event confidence

FastCall consumes only segment medians, one observation per segment. Segment
length, depth, IN/OUT composition and measurement error do not enter its
likelihood. Its truncated class supports have disjoint interiors. An interior
value therefore has posterior 1 for its sole supported class, irrespective of
support or noise. Shared endpoints are the exception and can have competing
class probabilities. Fragmentation changes priors and fitted widths without
changing those interior probabilities. Both numerical backends reproduce this.

`ProbCall` and VCF `FCP` are conditional class membership. VCF descriptions and
the manifest explicitly exclude calibrated event confidence. They must not be
used as empirical error rates. Numerical stability and coherent fitting do not
establish calibration.

An uncertainty-aware replacement needs a likelihood for noisy segment estimates,
including window dependence and purity uncertainty. Evaluate it using independent
truth-labelled samples, keeping patients and shared controls separate between
fit and validation sets. Measure reliability, Brier score and log loss by length,
depth, IN/OUT composition, noise and purity; include null segments and repeated
segmentations of the same event. Do not treat multiplied segment length as an
independent-observation count. No truth-labelled cohort was available for this
assessment, so empirical calibration remains open and no new confidence model
is introduced.

Regression evidence is in `tests/test_reporting_semantics.py`,
`tests/test_analysis.py` and `tests/test_signal_preservation.py`.
