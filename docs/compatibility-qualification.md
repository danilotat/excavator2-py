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
| Non-normal paired/pooling | Original fixtures and live CI, including BAM preparation | One complete fresh target-to-calls non-normal acceptance fixture |
| Target geometry | Eight cases; live old/new comparison; supplied target exact | Additional assemblies, unusual contigs and broader interval edge cases |
| Target features | Eight cases plus all supplied arrays; binary32/R decimal replay | UCSC 3,000-block algorithm-switch boundary, additional missing-data patterns |
| Read counting | Endpoints, overlaps, terminal/chunk cases; supplied exact counts | CRAM, empty selected streams and broader alignment-flag combinations |
| Numerical kernels | Original fixtures, native/reference comparison, sanitizer smoke checks | Broader tie/degenerate-arm cases and additional supported runtime configurations |
| Reproducibility | One/four-worker preparation tests; pinned oracle; installed wheels | CLI RNG-state interface for random FastCall ties; analysis parallelism remains unsupported |
| Plots and diagnostics | Scientific data tables/VCF comparisons | Plots remain unimplemented and unqualified; no visual parity claim |
| Packaging | Linux/macOS Python 3.11–3.13 CI and independent clean-wheel run | Release notices, migration guide and distribution qualification belong to M8 |

Next bounded step: **M6.2, a non-normal end-to-end fixture starting from target
configuration and reference files**, exercised through all three CLIs and compared
with the original in CI. Then address remaining numerical/feature boundary tests
and explicitly resolve plot scope before signing off M6. Keep scientific bug fixes
and performance changes separate from qualification.
