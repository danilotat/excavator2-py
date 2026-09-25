# M5 target generation

M5.3 is complete for the captured reference fixtures: Python geometry and
FASTA/BigWig features match the legacy outputs. The `target` CLI remains a
scaffold until artifacts and orchestration are integrated. M5.4 is next.

## Compatibility contract

The source of truth is `excavator2/lib/R/FilterTarget.R` at commit
`92f9c6e79fdeae32a09b99fa6994edc9a53e4efc`, executed unchanged in the pinned
legacy container. These are observed compatibility rules, including bugs;
correcting them belongs after call parity.

- Chromosomes are hard-coded as 1–22 and X. Prefix `chr` when the first BED
  chromosome name has more than three characters. Y and other contigs are ignored.
  Coordinate-file names are overwritten internally; numeric coordinate rows are
  interpreted in canonical chromosome order.
- On-target rows retain BED coordinates. They are not merged. Within a chromosome,
  BED order determines the gaps between chromosome start, target ends, target
  starts and chromosome end. Sorting the BED first would change the output.
- A gap is eligible when its length is at least `window + 400`. Flank is 200.
  With `n = floor(gap_length / window)`, boundaries are
  `gap_start + 200 + k * window`, for `k = 0..n`. Each OUT window starts one
  base after its left boundary and ends at its right boundary. Flanks are not
  subtracted before computing `n`: windows can overlap IN regions and extend
  beyond chromosome end.
- Chromosomes without BED entries instead use a grid from chromosome start to
  end, spaced by window size: OUT coordinates are `[grid[k], grid[k+1]-1]`.
  This differs from the populated-chromosome convention.
- IN rows are concatenated before OUT rows, then sorted numerically by start
  within canonical chromosome order. The captured duplicate-start IN rows retain
  their original order.
- Gap input skips its first line and uses columns 2–4 for chromosome/start/end.
  For bare chromosome names, the first `chr` occurrence is removed from gap names.
  Rows with `_` in their chromosome are excluded.
- A target row is removed if either endpoint lies in `[gap_start, gap_end)`.
  Full interval overlap is not tested: a row spanning a gap can survive, an end
  exactly at gap start is removed, and a start exactly at gap end survives.
- IDs `a1..aN` are assigned after filtering. The `MyTarget` character matrix has
  columns chromosome, start, end, ID, IN/OUT. Its values equal `Filtered.txt`.

The original wrapper also enforces a minimum window size of 10, checks reference
paths and only the first BED row's start/end, and nests output under
`<output>/<assembly>/<target>/w_<window>`. CLI/config translation is deferred.

## Captured evidence

All probes use window 100 and chromosome bounds 0–3000. Inputs, full output
arrays, stderr and SHA-256 provenance are in `tests/fixtures/legacy-target/`.

| Case | Original result | Behavior covered |
| --- | --- | --- |
| standard | 666 rows | IN/OUT coordinates, ignored Y, flanks, chromosome bounds |
| bare_names | 666 rows | Equivalent geometry without `chr` prefix |
| unsorted_overlap | 678 rows | Input order, overlapping regions, duplicate starts |
| gap_endpoints | 662 rows | Half-open endpoint filtering and spanning intervals |
| coordinate_names_ignored | 666 rows | Reversed coordinate labels leave results unchanged |
| no_eligible_gap | Fails | Empty eligible-gap loop produces nonfinite sequence bound |
| no_alternate_gap | Fails | Negative empty index drops all gap rows |
| missing_chromosome_gap | Fails | Empty chromosome gap loop indexes beyond bounds |

Expected failures are evidence, not desired new behavior. M5.2 should preserve
rejection of these cases with explicit Python errors; matching R error text is
not required. Broader edge-case coverage will be added with implementation.

Fixture tests run in the existing Linux/macOS wheel CI matrix. They compare every
Python geometry output cell, row order and chromosome list against the captured
original outputs, and check failure behavior. This is fixture-based geometry
concordance; live feature comparison is now also enabled. Full-target
concordance remains M5.4. Existing live
preparation/analysis concordance remains enabled.

The internal entry point is `excavator2.target.target_geometry(bed, coordinates,
gaps, window)`, returning a five-column NumPy string matrix. It creates no files
and performs no reference feature extraction. Geometry and compatibility policy
remain Python; NumPy handles array generation, sorting and filtering. Input scope
is tab-separated files with integer genomic coordinates, a nonempty BED and at
least 23 canonical coordinate rows. Extra BED columns are ignored. Full R text
parser equivalence, arbitrary assemblies and performance are not qualified.

## Remaining bounded steps

1. **M5.2 — Python geometry (complete):** readable Python/NumPy implementation;
   exact comparisons against all five success fixtures and explicit errors for
   all three failure cases. Separate from reference extraction.
2. **M5.3 — reference feature contract and extraction (complete for fixtures):** characterize the original
   FASTA/BigWig path with tiny references, including missing/ambiguous bases,
   uncovered values, boundaries, and decimal serialization. GC/MAP extraction
   uses BED `[target_start-1, target_end)`; first-base extraction uses
   `[target_start, target_start+1)`. The original map parser saves column 6 of
   `bigWigAverageOverBed` output. Prove feature parity before integration. Prefer
   existing native-backed readers; add custom C++ only for a measured bottleneck.
3. **M5.4 — integration and acceptance:** wire YAML/config, target artifacts and
   the `target` CLI into preparation; add live target comparisons to CI and run
   all-new `target → prepare → analyze` on the supplied data against original
   calls. M5 is complete only after this acceptance run passes.

M5.3 does not qualify arbitrary assemblies, full target CLI behavior, performance,
or end-to-end calls from newly generated targets. Reference features are qualified
on the fixtures described below.

## M5.3 reference features

`excavator2.target_features.target_features(target, fasta, bigwig)` now returns
per-chromosome `gc`, `mappability` and `first_base` arrays. FASTA access uses pysam;
BigWig access uses pyBigWig, which was already a package dependency. Policy and
rounding remain readable Python; reference I/O uses those native-backed libraries.
No custom C++ or CLI integration was added in this step.

The implementation uses `stats(..., type="mean", exact=True)`: uncovered bases
are excluded from the mean, and approximate zoom summaries are bypassed, as
specified by the [pyBigWig documentation](https://github.com/deeptools/pyBigWig#compute-summary-information-on-a-range).
The unchanged legacy scripts confirm these additional rules:

- MAP uses covered-base mean (column six), with zero for no coverage, including
  chromosomes with no BigWig data. Values pass through six significant decimal
  digits, reproducing `bigWigAverageOverBed` text output before R reads it.
- GC counts G/C case-insensitively over the full requested sequence length,
  including ambiguous bases in the denominator. Values pass through six decimal
  places, reproducing `bedtools nuc` text output.
- FRB preserves base case and ambiguous characters. It returns two string columns:
  target start and the base at that zero-based offset.
- The shell selects chromosomes using `grep -w` on entire target rows. Bare names
  such as `1` can therefore select rows from other chromosomes whose coordinate
  columns contain `1`. This bug remains explicit in the Python implementation.
- GC intervals extending past the FASTA contig end are skipped; BigWig means can
  still be computed from covered bases within bounds. FRB is skipped when its
  one-base interval is outside the FASTA. If no FRB rows survive for a chromosome,
  the original fails and Python raises an explicit error. Feature array lengths
  can consequently differ. No silent padding or truncation is performed.
- A zero target start produces a negative GC/MAP BED coordinate and fails.

Seven cases in `tests/fixtures/legacy-target-features` cover fractional GC and
mapping values, lowercase and ambiguous bases, partial/missing/zero BigWig
coverage, bare chromosome names, and reference boundaries. Three cases succeed
with exact array comparisons; four reproduce legacy failures. Small FASTA,
index, BigWig, target inputs, expected arrays, logs and hashes are committed.

Reproduce live original-versus-current feature comparison with:

```sh
python tools/oracle/characterize_features.py \
  --output .oracle/feature-concordance --compare-current
```

Use a fresh output directory with Docker available. The harness extracts the
pinned source commit and runs its unchanged `TargetCreate.sh` and R save scripts.
The existing CI concordance job now performs this comparison using the installed
current wheel and retains the results alongside preparation/analysis reports.
Local feature parity is established on these fixtures; full supplied-reference
qualification and `target → prepare → analyze` acceptance remain M5.4. Arbitrary
FASTA/BigWig variants, shell metacharacters in names and performance are not yet
qualified.
