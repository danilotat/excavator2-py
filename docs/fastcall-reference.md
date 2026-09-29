# FastCall: Python control and scalar C++ kernels

FastCall consumes one log2-ratio value per segment. Python owns initialization,
iteration, cellularity correction and labels; the repeated numerical calculations
have both NumPy/SciPy and scalar C++ implementations.

```python
from excavator2.fastcall import correct_cellularity, fit_fastcall, assign_labels

values = correct_cellularity([0, -0.4, 0.4, 1, 2, 3], cellularity=1)
fit = fit_fastcall(values)  # backend="python" selects the readable implementation
calls = assign_labels(fit.posterior)
```

Preserve the original values separately for copy-number reporting. See
[M3 analysis](hslm-analysis.md) for CLI integration.

## Model and fitting

The five fixed means are `[-3, -1, 0, 0.58, 1]`. Class supports are inclusive;
`upper` and the positive magnitude `lower` define the normal interval
`[-lower, upper]`. Fitting and reporting use the same truncated densities.

The M-step optimizes the truncated-normal objective, with deviations at least
`0.001`. A deviation of positive infinity denotes the uniform limit within that
class's bounds. Empty components receive zero prior and retain their deviation.
The E-step normalizes log densities; unsupported observations raise instead of
receiving a nearest-mean guess. See [A19](fastcall-a19.md) and
[A20/A21 and the convergence correction](fastcall-fitting.md) for the rationale.

`fit.trace` contains five deviations, five priors and the actual truncated
mixture log-likelihood after each M-step. Convergence requires stable likelihood,
precisions and priors. `fit.converged` is false if the 1,000-iteration cap is
reached. Analysis manifests identify the posterior as `truncated` and the trace
statistic as `truncated_log_likelihood`.

Inputs must be finite, nonempty 1-D values within the outer bounds `[-20, 20]`.
Thresholds must satisfy `0 < upper < 0.9`, `0 < lower < 1.3`; cellularity must be
positive and no greater than one. A class probability of 1 inside its exclusive
support expresses the hard boundary assumption, not calibrated biological certainty.

## Numerical implementation

- `src/excavator2/fastcall.py`: fit orchestration and label assignment.
- `src/excavator2/reference/fastcall.py`: readable E/M steps and log probabilities.
- `cpp/fastcall.cpp`: the corresponding scalar kernels.
- `cpp/bindings/module.cpp`: validation, allocation and GIL release.

The public fit API converts input to aligned contiguous binary64 arrays. The
private native boundary rejects implicit copies, invalid shapes and parameters.
Vectors have shapes `(n,)` or `(5,)`, probabilities `(n, 5)` and bounds `(5, 2)`;
matrices are row-major. E-step returns probabilities and log-likelihood; M-step
also accepts bounds. Positive-infinite deviations are accepted only by the
truncated E/M kernels. Inputs are not mutated; outputs own their storage.
Do not mutate shared inputs during native calls.

## Label ties

Label assignment retains the R `max.col` relative tie tolerance of `1e-5`.
Without an RNG state, near ties raise explicitly. For replay, pass the original
626 signed integers of R's Mersenne-Twister `.Random.seed` to `assign_labels`.
The returned `calls.r_seed` continues that stream. A scalar seed or NumPy seed
is not an equivalent state. See the [R max.col documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/maxCol.html).

## Verification

The tests cover independent truncated-density and numerical-integration oracles,
objective improvement and stationarity, both variance limits, empty components,
extreme tails, boundary support, custom thresholds, convergence, input validation
and concurrency. Original fixtures remain useful inputs and historical evidence;
the corrected fit is not required to reproduce their old parameters or calls.

The standalone kernel smoke test passes AddressSanitizer and UndefinedBehaviorSanitizer;
see `cpp/tests/README.md`. Run `uv run pytest` for the full suite.

A local five-repeat benchmark after A20/A21 used 10,003 segments on macOS ARM64,
Python 3.12.6 and NumPy 2.5.3: median full-fit time was 9.27 ms for Python and
0.95 ms for native (9.7x), with two iterations. This excludes HSLM, I/O and labels
and is not a full-pipeline speedup claim. Reproduce with
`python benchmarks/fastcall.py --segments 10003 --repeats 5`.
