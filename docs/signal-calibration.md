# Signal preservation before calling

P14 ([#14](https://github.com/danilotat/excavator2-py/issues/14)), P13
([#13](https://github.com/danilotat/excavator2-py/issues/13)) and A17
([#5](https://github.com/danilotat/excavator2-py/issues/5)) change the preparation
and ratio contract. Legacy normalized artifacts cannot recover erased signal and
are rejected by analysis. Re-run preparation from BAMs. There is no legacy mode.

## Supported domain

Analysis requires matched windows, nonnegative integer observed counts, a diploid
reference, independently established relative sample exposures, and a shared
multiplicative technical response to GC, mappability and length across samples.
Exposure means expected reads per copy at a fixed window before its shared
technical response. Use external standards or independently validated neutral
calibration regions; do not estimate exposure from the median of the analyzed
IN/OUT windows or total reads when genome-wide CN can change. All exposures must
use the same scale. A common multiplicative change to all exposures cancels.

The user supplies these scientific assertions; the software validates their
presence and values but cannot establish their truth from read counts. Unknown
baseline/ploidy, unmatched sample-specific covariate bias, and profile-derived
exposures are unsupported. These are identifiability limits, not effects that a
more robust per-sample median can always solve. A shared response is deliberately
cancelled by comparison rather than estimated separately in each sample.

```yaml
# calibration.yaml (names must exactly cover the selected tests and controls)
baseline: diploid-reference
bias: shared
exposure_source: independent
exposures:
  Control1: 1.0
  Test1: 1.2
```

The numbers above are illustrative, not calibration estimates for a dataset.

```sh
excavator2 prepare --samples bams.yaml --target target/ --output prepared/
excavator2 analyze --samples design.yaml --input prepared/ --target target/ \
  --output results/ --experiment paired --calibration calibration.yaml
```

## P14: zero observations and missing data

Preparation preserves zeros, including an all-zero class. It never uses another
window's minimum. Analysis adds a fixed 0.5 **read** to every raw window count
before dividing by exposure. For test counts t and controls c_j with exposures
s_t and s_j, the ratio is:

`log2(((t + 0.5) / s_t) / mean_j((c_j + 0.5) / s_j))`.

Matched window lengths cancel. Pooling gives equal weight to controls after
exposure calibration. Computation is in log space to avoid overflow in scaling.
The half-read prior is depth-aware: a test zero against 100 reads is a larger
loss than a test zero against one read. It provides finite regularization, not
proof of a biological deletion, exact zero copy number, or a depth uncertainty
model. At low depth it attenuates ratios, including cancellation of technical
bias; interpretation requires adequate coverage. Changing an unrelated low-depth
window cannot change this ratio.

Any zero in any selected control is unsupported, including both-zero windows and
all-zero reference classes. Analysis fails the whole run before publishing
results; it does not silently discard those windows or substitute a pooled
positive control. Missing/non-finite/negative/fractional counts are rejected,
never recoded as observed zero. A zero must represent a measured window, not a
missing assay or failed extraction. The existing BAM counting/filtering defects
listed in [preparation](preparation.md) are separate and unchanged by this fix.

## P13: covariate-associated signal

Preparation now records uncorrected density (`counts / (end - start)`) for the
seven-column matrix and retains raw integer counts as the analysis input. It
performs no per-sample size/MAP/GC bin rescaling. A biological gain confined to a
GC, MAP or length bin therefore survives. Shared technical effects cancel in the
matched test/reference ratio, up to the explicit half-read regularization.
Sample-specific technical effects remain confounded with biological effects and
are outside the supported domain. Legacy correction checkpoints are no longer
emitted; the manifest identifies `depth_policy: raw-counts-v1`.

## A17: baseline and broad alterations

There is no median centering, separately by class or jointly. A majority-altered
IN class with neutral OUT, or the reverse, retains its calibrated gain and leaves
neutral windows neutral. Whole-profile biological scaling also remains in the
ratios when independent exposures support that interpretation; a corresponding
change in sampling exposure is removed by the supplied exposure factors.
Without independent calibration those explanations cannot be distinguished, so
analysis requires calibration rather than inferring absolute CN silently.

This contract preserves signal **before** HSLM and FastCall. A perfectly constant
profile still fails HSLM's positive-variance requirement; finite ratios do not
promise a supported downstream fit or clinical validity. See [copy-number interpretation](copy-number-interpretation.md) for corrected
cellularity reporting, diploid-equivalent semantics, pool weighting and the
remaining probability-calibration limitation.

The analysis manifest records the full calibration and its file hash, raw depth
policy, half-read pseudocount and `centering: none`. Existing calibrated ratio
checkpoints remain available. Regression tests in `test_signal_preservation.py`
measure retained biological effect and residual shared technical bias separately,
and cover zero/reference failures, missing values, class-majority shifts, whole
profile shifts and exposure-scaled pooling. Historical normalization and
end-to-end legacy concordance claims no longer apply to this path; counting and
isolated segmentation fixtures remain useful.

CI now runs the calibrated preparation/analysis regression suite against the
installed wheel and retains live target and isolated HSLM/arm evidence checks.
Legacy whole-pipeline equality gates and their automatic benchmark publication
were retired because normalization and baseline semantics intentionally differ.
The historical oracle scripts/reports remain available for reproducing the old
contract; they are not acceptance gates for the calibrated path.
