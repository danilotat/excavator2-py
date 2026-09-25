# M6 compatibility qualification

M6.1 checks the M5 implementation at commit
`ac53247279f28966437c7788c49caa0196064941`. This is a bounded qualification step,
not completion of the entire compatibility milestone or permission to repair
legacy scientific behavior.

## Remote evidence

[Python port run 36102519608](https://github.com/danilotat/excavator2-py/actions/runs/36102519608)
passed all seven jobs for that exact commit:

- Linux (`ubuntu-latest`), Python 3.11, 3.12 and 3.13.
- macOS (`macos-latest`), Python 3.11, 3.12 and 3.13.
- Live original-versus-current concordance for geometry, reference features,
  BAM preparation, paired analysis and pooled analysis.

Each package job runs lint, formatting, fixture tests, native sanitizer smoke
checks, wheel builds and tests against the installed wheel outside the checkout.
The live job runs the digest-pinned original runtime and current wheel on fresh
small inputs. Success establishes the configured runner matrix; it does not imply
Windows, every CPU architecture, every dependency version or alternate SIMD paths
are qualified.

## Independent clean-wheel acceptance

A new Python 3.12 environment was populated from a freshly built wheel and its
runtime dependencies, with no editable installation. Package and extension paths
were verified under that environment's `site-packages`. Commands ran from outside
the checkout with Python `-I` isolation. Only input paths in the example YAML were
resolved to absolute paths; scientific settings and input files were unchanged.

The verification covers all 137 tests and a new supplied-data
`target → prepare → analyze` paired run. Target arrays compare directly with the
converted original target; preparation arrays compare with the M4 artifacts
already qualified against the original; final analysis compares directly with the
saved original results. Reports and runtime/build provenance are retained in
`tools/oracle/m6-clean-wheel-report.json`.

To reproduce, create a fresh environment, install the built wheel plus pytest,
resolve input paths for execution outside the checkout, and run:

```sh
/path/to/clean/bin/python -I -m pytest /absolute/repository/tests
/path/to/clean/bin/python -I -m excavator2 target --config config.yaml --output target/
/path/to/clean/bin/python -I -m excavator2 prepare --samples samples.yaml --target target/ --output prepared/ --threads 4
/path/to/clean/bin/python -I -m excavator2 analyze --samples design.yaml --target target/ --input prepared/ --output results/ --experiment paired
```

Use `tools/oracle/compare_target.py` and `compare_analysis.py` for the exact target
and decision-gated analysis comparisons. Compare every preparation NPZ array,
not just final calls. See the M5 reports for the original comparison basis.

## Qualification matrix and remaining work

| Area | Evidence established | Remaining gap |
| --- | --- | --- |
| Supplied paired pipeline | New target through final calls; 281,567 windows per sample | Supplied example has no non-normal calls |
| Non-normal paired/pooling | Original fixtures, stage CI and M6.2 complete fresh-reference chain | Remote confirmation of the newly added full-chain CI step |
| Target geometry | Eight cases; live old/new comparison; supplied target exact | Additional assemblies, unusual contigs and broader interval edge cases |
| Target features | Eight cases plus all supplied arrays; binary32/R decimal replay | UCSC 3,000-block algorithm-switch boundary, additional missing-data patterns |
| Read counting | Endpoints, overlaps, terminal/chunk cases; supplied exact counts | CRAM, empty selected streams and broader alignment-flag combinations |
| Numerical kernels | Original fixtures, native/reference comparison, sanitizer smoke checks | Broader tie/degenerate-arm cases and additional supported runtime configurations |
| Reproducibility | One/four-worker preparation tests; pinned oracle; installed wheels | CLI RNG-state interface for random FastCall ties; analysis parallelism remains unsupported |
| Plots and diagnostics | Scientific data tables/VCF comparisons | Plots remain unimplemented and unqualified; no visual parity claim |
| Packaging | Linux/macOS Python 3.11–3.13 CI and independent clean-wheel run | Release notices, migration guide and distribution qualification belong to M8 |

Next bounded step: **M6.3, target-feature boundary qualification**, starting with
the UCSC 3,000-block algorithm switch and its rounding implications. Then address
remaining numerical edge cases and explicitly resolve plot scope before signing
off M6. Keep scientific bug fixes
and performance changes separate from qualification.


## M6.2 fresh-reference non-normal chain

The installed-wheel run now covers the complete original and current target,
preparation and analysis CLIs on generated inputs: 23 reference contigs, 4,113
windows, two controls and one test with decreases/increases on chr1 and chr2.
Both paired and pooling produce four non-normal calls (losses and gains).

All target/GC/MAP/FRB arrays, all 12,339 raw counts and all prepared matrices match
exactly. Each mode's four scientific output files matches exactly apart from VCF
fileDate; this fixture has zero measured continuous HSLM drift. The harness rejects
missing output files, fewer/different calls, and loss of either loss/gain coverage.
It requires exact target and preparation equality before accepting final calls.

See `tools/oracle/m6-end-to-end-report.json` and reproduce using:

```sh
/path/to/installed-wheel/python tools/oracle/end_to_end.py --output .oracle/end-to-end
```

The current path generates its own target; converted original artifacts are used
only as comparison inputs. No downloaded genome/BAM assets or previous oracle
runs are required. CI now invokes the same harness and retains its complete logs,
inputs and report. Local acceptance passed; the new remote CI step is pending.
Six comparator tests additionally prove that valid-checksum scientific changes
in target features, row ordering, reference bases and centromeres are rejected.
The full local suite has 143 passing tests.

Plot qualification remains open: the original plotting code cannot obtain UCSC
chromInfo for the synthetic assembly in the offline container and fails in both
modes. This is explicit in the report; the pass applies to scientific data and
calls only. A first fixture with longer overlapping IN windows also exposed
nonmonotonic midpoints and failed in both original inference and port validation.
The passing fixture uses shorter intervals, without changing either scientific
implementation. Reproducer details are retained in the local improvement log.
