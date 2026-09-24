# Python FastCall reference

The first numerical port is available through the Python API. The command-line
pipeline still stops explicitly: HSLM, segment-table construction, preparation,
and output writers are not integrated yet. No C++ numerical kernel was added in
this step.

## Readable implementation

- `src/excavator2/fastcall.py`: initialization, cellularity correction, iteration
  and stopping policy, assignment, and optional replay of the original R RNG state.
- `src/excavator2/reference/fastcall.py`: normal densities, truncated E-step,
  fixed-mean M-step, posterior normalization, and underflow fallback. These small
  NumPy/SciPy functions are the references for future C++ kernels.

The API consumes **one log2-ratio value per segment**, not the full window vector.
Callers must preserve the original uncorrected values for later CN reporting.
For example, from the repository root after installing the development package:

```python
import numpy as np
from excavator2.fastcall import correct_cellularity, fit_fastcall, assign_labels

with np.load("tests/fixtures/legacy-fastcall/expected.npz", allow_pickle=False) as oracle:
    original = oracle["mdata"]
    calling_values = correct_cellularity(original, cellularity=1.0)
    fit = fit_fastcall(calling_values, upper=0.35, lower=0.5)
    calls = assign_labels(fit.posterior, r_seed=oracle["seed_before"])
    print(calls.labels)
    print(fit.iterations)
```

`fit.trace` records five deviations, five priors, and the legacy stopping statistic
at each iteration. Means remain fixed. The initial and E-step bounds differ, as
in the original. Empty-component priors, inclusive truncation endpoints, tiny-SD
retention, underflow fallback, and the 1,000-iteration cap are preserved.

The stopping statistic reproduces R's column-major vector recycling in
`sum(PosteriorP(...) * prior)`. It is deliberately not a conventional mixture
likelihood. Fixing that expression would change the algorithm and is deferred.

## Random ties

R's default `max.col` treats sufficiently close values as ties using a relative
1e-5 threshold and randomized assignment. See the [R documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/maxCol.html)
and [R 4.4 implementation](https://github.com/wch/r-source/blob/R-4-4-branch/src/appl/maxcol.c).
The Python reference implements the same sequential reservoir decision policy.

For exact replay, pass the 626-element `.Random.seed` exported from R's
Mersenne-Twister generator. The NumPy MT engine is initialized with those state
words; each R uniform consumes one raw 32-bit word. Ordinary NumPy random floats
or the same integer seed do not establish R compatibility. The returned
`calls.r_seed` can resume the stream for a later assignment batch. Other R RNG
kinds are rejected. See the [R RNG state implementation](https://github.com/wch/r-source/blob/R-4-4-branch/src/main/RNG.c).

Without an R state, unique maxima can be called, but a near tie raises a clear
error. This avoids silently substituting a first-maximum policy or inventing a
historical seed. That path does not advance or return an R RNG state. Intermediate
ties that later lose still consume random draws when replay is requested.

## Verification and limits

The default test suite compares against the original pinned R code for:

- The original five-state fixture.
- Unequal state populations, with a segment count not divisible by five.
- Truncation endpoints, including a 40-iteration case.
- All-normal, single-segment, extreme-value, and nondefault-bound cases.
- Segment values reconstructed from the saved real paired HSLM result.
- Exact and near ties, RNG-state continuation, and cellularity correction.

Labels, iteration counts, and final RNG states match exactly on these fixtures.
On the local Python 3.12/macOS ARM64 run, the largest absolute iteration-trace
error was approximately 7.1e-15; the largest posterior error was approximately
2.3e-16. Parameter/trace tests use rtol=1e-13 and atol=1e-14; posterior/probability
tests use rtol=1e-13 and atol=2e-15. These tolerances cover observed floating-point
drift and never allow a different discrete decision or stopping iteration.
Cross-platform CI still needs to qualify these bounds independently.

The trace fixture generator verifies that instrumented and uninstrumented R fits
are identical. Existing goldens were retained; new cases are in
`tests/fixtures/legacy-fastcall/edges/`, with script/source hashes and binary values.
`tools/oracle/characterize_fastcall_edges.R` and `tools/oracle/pack_fastcall.py`
provide regeneration tools. The packer refuses to overwrite an existing fixture
directory and rejects unexpected missing values.

The current fit API requires finite, nonempty 1-D data, positive cellularity up to
one, and ordered model bounds (`0 < upper < 0.9`, `0 < lower < 1.3`). It fails
explicitly on non-finite iteration results; it does not add numerical stabilization
to repair legacy failures. This is fixture-backed numerical parity for FastCall,
not proof of universal equivalence or an end-to-end CNV calling release.

Next: add the focused C++ E/M and posterior kernels behind this Python control
logic, requiring the same fixtures and decision checks to pass.
