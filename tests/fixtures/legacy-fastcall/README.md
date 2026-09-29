# Original FastCall: five-state fixture

`expected.npz` contains binary arrays exported from the original `LibraryFastCall.R`
in the pinned image. Load with `numpy.load(..., allow_pickle=False)`. It includes
input segment values, bounds, fitted parameters, iteration count, posterior matrix,
labels/probabilities, and R RNG states before and after the call.

`calls.tsv` is a readable rendering; use the NPZ for numerical comparisons.
`manifest.json` records provenance, script hashes, seed, dtypes, shapes, and archive
checksum. The fixed seed belongs only to this synthetic characterization; the full
BAM-based baseline runs do not inject a seed.

Regenerate only deliberately using `tools/oracle/characterize_fastcall.R` against
the frozen source and `tools/oracle/export.R`. Preserve binary64 values when converting
column-major R exports to NumPy. Do not overwrite expected values to silence a
future differential-test failure. This fixture exercises all five states but does
not replace the real acceptance dataset or establish biological accuracy.


Additional fixtures in `edges/` cover unequal populations, truncation endpoints,
empty components, a single segment, extreme values, nondefault bounds, the real
baseline segmentation, cellularity, and random near ties. They include per-iteration
traces. The R generator checks identical fits with tracing enabled and disabled.
`edges/manifest.json` records their provenance and archive hashes.
