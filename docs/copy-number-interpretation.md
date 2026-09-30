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

The calibration assertion `diploid-reference` applies to **every analyzed
window in every control**, including sex chromosomes and PAR. The contaminating
component in the cellularity inversion must also match that baseline. Chromosome
names cannot establish these facts. Absolute sample ploidy, sex, PAR boundaries
and allele counts are not inferred or accepted as metadata by this pipeline.
Non-diploid or unknown reference assertions are rejected. A falsely declared
diploid reference cannot be detected from these counts alone.

| Reference domain | Interpretation / support |
| --- | --- |
| Validated diploid autosome | A twofold ratio gives DECNF 4 at p = 1 |
| Haploid X outside PAR | Unsupported reference; a twofold ratio would give diploid-equivalent 4, not absolute 2 |
| Validated two-copy PAR in all controls | Satisfies the diploid baseline assertion; no automatic PAR detection |
| PAR with unknown or mixed reference dosage | Unsupported; chromosome-level sex alone is insufficient |
| Three-copy aneuploid reference | Unsupported; a twofold ratio would give equivalent 4, not absolute 6 |

This addresses misleading absolute naming, not general ploidy-aware inference.
Adding absolute CN later requires region-specific reference and contaminant
ploidy, sample context, validated PAR handling, and an appropriate calling model.
The fixed FastCall means and class thresholds retain their diploid-relative meaning.

## P15 assessment: equal weighting is already implemented

The current `ratios` function averages `(control_count + 0.5) / exposure` equally
across controls. Independent exposure calibration removes the original depth
weighting. In the audit-shaped probe (100 reads, one flat control, another with
a twofold gain in ten windows), the affected ratio is `log2(100.5/150.5)`.
Increasing only the second control's depth and exposure tenfold gives
`log2(100.5/150.275)`. The difference is below 0.004 log2 units across all windows;
it comes from applying the half-read prior before exposure scaling. Exact scale
invariance is intentionally not claimed at low depth.

The control gain still distorts the reference: equal weighting does not remove
control CNVs. There is no quality weighting, outlier rejection, or robust
aggregation. Controls must satisfy the stated baseline; automatic CNV handling
and alternative weights require validation on independently characterized
controls. No pooling algorithm change is justified by this probe alone.

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
