# A20/A21: coherent fitting and stable probabilities

This change addresses [A20 (#8)](https://github.com/danilotat/excavator2-py/issues/8)
and [A21 (#9)](https://github.com/danilotat/excavator2-py/issues/9). It also replaces
the invalid stopping statistic described in [A08 (#4)](https://github.com/danilotat/excavator2-py/issues/4)
so convergence is measured against the objective being fitted. There is no legacy mode.

## Width update

For one component with fixed mean mu, responsibility sum W and weighted mean
squared deviation m2, the objective per unit weight is:

```text
Q(sigma)/W = -log(sigma) - log(sqrt(2*pi)) - m2/(2*sigma^2) - log Z(sigma)
Z(sigma) = Phi((upper-mu)/sigma) - Phi((lower-mu)/sigma)
```

The former RMS update omitted the derivative of `log Z`. The corrected update
solves `E[(X-mu)^2 | lower <= X <= upper] = m2` by bisection in log sigma.
In precision `t = 1/(2*sigma^2)`, the derivative is `E[(X-mu)^2] - m2` and
its derivative is minus the variance of `(X-mu)^2`. The objective is concave
in t, so the stationary point is its maximum.

The constrained lower width is `0.001`, matching the initialization floor and
preventing variance collapse on identical observations. At the other limit,
if m2 is at least the uniform second moment, the optimum is the uniform density
on the finite support. This is represented by `sigma = +inf`, with no arbitrary
large-width cutoff. Empty components retain their width and receive zero prior.
Candidate updates are accepted only when their evaluated component objective
does not decrease. These changes are shared by both numerical backends.

For A20's `[0.2, 0.3]` example, mean 0 and bounds `[-0.5, 0.35]`:

| Update | Width | Summed truncated log density |
| --- | ---: | ---: |
| Starting point | 0.4 | 0.29164932832 |
| Original RMS update | 0.25495097568 | 0.12818861934 |
| Corrected update | 1.59228318 | 0.32520182361 |

## Probability calculation

Responsibilities use log weights and normalization after subtracting the largest
log weight. Out-of-support components have log weight `-inf`; zero priors remain
zero. Rows with no supported positive-weight component fail explicitly.

Tail normalization uses scaled Gaussian tails, avoiding subtraction of CDF values
that have both rounded to one. Python uses [SciPy erfcx](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erfcx.html);
C++ uses erfc for moderate arguments and the Mills-ratio asymptotic expansion
for large arguments. Densities are expressed relative to the largest density
inside the interval, avoiding subtraction of enormous negative log densities.
Narrow standardized intervals and near-uniform second moments use 16-point
Gauss-Legendre integration where their integrands are smooth.

The `[0.30]` probe with all widths `0.001` now assigns normal probability 1,
consistent with its support. In the independent untruncated `[0.29]` probe with
widths `[0.01, 0.01, 0.001, 0.001, 0.01]`, amplification wins as the log-density
calculation requires. No nearest-mean or `Inf -> 100` fallback remains.

## Convergence and diagnostics

The E-step returns the mixture log-likelihood as well as responsibilities.
After each M-step, the likelihood must not decrease beyond relative numerical
tolerance `1e-10`. Stopping requires likelihood change within `1e-8*(1+abs(previous))`,
relative precision change within `1e-8`, and prior change within `1e-8`.
Precision represents the uniform limit as zero. The 1,000-iteration cap remains;
reaching it returns `converged=False`.

The last column of `fit.trace` is now truncated log-likelihood, not the previous
recycled posterior/prior statistic. Analysis manifests record this change.
Positive infinity in a deviation column is a valid uniform component, not a
numerical error. Calls and fitted parameters can change from the A19-only version.

## Verification

Both backends are checked against [SciPy truncated-normal densities](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.truncnorm.html),
finite-difference objective derivatives and independently scaled adaptive
integration. Coverage includes empty components, uniform and floor optima,
positive/negative tails, widths from `0.001` to `1e8`, tiny priors, exact boundaries,
unsupported observations and permutation invariance. Counts of 5, 10 and 15
exercise the former one-iteration stopping failure. Native buffers also pass
AddressSanitizer and UndefinedBehaviorSanitizer.

These fixes do not calibrate biological uncertainty. Away from a shared endpoint,
the hard component supports still allow only one class and therefore produce
probability 1.

Compared with the A19-only commit `aa450ad`, the six supported saved FastCall
fit cases change one label: in `nondefault`, the shared threshold `+0.2` changes
from normal to duplication (winning probability about `0.587402`). The other
two fit cases remain unsupported because they include values outside `[-20, 20]`.
All ten saved analysis profiles and the supplied sample retain their labels.
A full rerun from the supplied prepared artifacts reproduces its 281,567 HSLM
rows and FastCall/VCF output text, excluding the VCF date. This sample has no
non-normal calls. Local comparison reports are in `.oracle/a20-a21/`.
