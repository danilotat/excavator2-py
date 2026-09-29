# A19: final probabilities respect fitted class support

On `dev/improvements`, FastCall uses the E-step's truncated densities to compute
final probabilities. `main` is the unchanged legacy baseline. This correction
addresses [A19](https://github.com/danilotat/excavator2-py/issues/7).

## Behavior

Previously, the E-step restricted each component to its interval, but reporting
discarded those restrictions. A broad amplification curve could win even for a
negative value; another component could also steal a value inside amplification
support. The final call could therefore contradict the model used for fitting.

After the last M-step, FastCall recomputes probabilities with the same
E-step kernel, final fitted parameters, and inclusive bounds. It does not reuse
responsibilities from before the last parameter update. Both Python and C++
backends use their existing kernels. No original R code or golden fixture changes.

With the default thresholds, supports are `[-20, -1.3]`, `[-1.3, -0.5]`,
`[-0.5, 0.35]`, `[0.35, 0.9]`, and `[0.9, 20]`. Custom `d`/`u` thresholds apply
to final calls too. Exactly shared endpoints allow both adjacent components;
immediately outside an endpoint a component has zero probability.

For a class j, its reporting weight is now:

```text
prior[j] * normal_pdf(x; mean[j], sd[j]) * I(lower[j] <= x <= upper[j])
----------------------------------------------------------------------------
 normal_cdf(upper[j]; mean[j], sd[j]) - normal_cdf(lower[j]; mean[j], sd[j])
```

The five weights are normalized to sum to one. Outside shared endpoints only
one class is allowed, so its probability is 1. This is a consequence of the hard
support assumption, not calibrated certainty about biological copy number.

## Scope

The API and CLI use the correction directly, without a compatibility option.
The analysis manifest records `fastcall_posterior: truncated`.

Initialization, M-step, stopping statistic, iteration count, fitted parameters,
and trace are unchanged. A20's truncated M-step objective and A08's stopping
statistic remain open. A21's ordinary-space arithmetic is also unchanged: the
corrected reporting path rejects non-finite probabilities and any fallback that
places mass outside support. This includes values beyond the outer bounds and
the `[0.3, 0.3]` underflow example. Analysis fails without publishing partial results.
Fallbacks inside support and probability calibration are not solved by this patch.

## Evidence

The regression tests run both backends and retain the legacy goldens:

| Synthetic input | Value | Legacy label / probability | Corrected label / probability |
| --- | ---: | --- | --- |
| Ten zeros, -0.4, 0.4, 1, 2, 3, 4, 5 | -0.4 | Amplification / 0.823040679649 | Normal / 1 |
| Fifty repetitions of [-0.3, 0.3], then 0.91, 5 | 0.91 | Normal / 0.825767131 | Amplification / 1 |

Each example changes exactly one label. Other probabilities can change. These
examples test class consistency under the configured thresholds, not biological
ground truth or the population frequency of errors.

The comparison on 2026-09-29 covered all eight saved FastCall fit cases, ten
saved analysis profiles (paired, pooling and chromosome-arm cases), and the
supplied sample. Six FastCall cases and all eleven analysis/sample profiles
retained their labels. Probability vectors changed by more than 1e-12 for 28
segments across those six FastCall cases and eight across the arm profiles.
The two remaining FastCall cases (`boundaries` and `extreme`) contain values
outside `[-20, 20]` and are explicitly rejected. The original fixtures remain
as historical evidence.

A full analysis rerun from the supplied prepared artifacts retained all 281,567
HSLM rows, its normal segment, and the FastCall/VCF output text (excluding VCF
date). This sample has no non-normal calls, so the synthetic regressions are
necessary to demonstrate the correction. Local detailed comparisons are saved
under `.oracle/a19-correction/`.

Boundary tests compare against independent SciPy truncated-normal densities,
including custom thresholds and neighboring representable floating-point values.
API and CLI tests cover corrected calls, output probabilities, and unsupported
numerical fallback rejection. RNG replay is tested separately.

Run the focused checks with:

```sh
uv run pytest tests/test_fastcall.py tests/test_fastcall_native.py tests/test_analysis.py
```
