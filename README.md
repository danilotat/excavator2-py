# EXCAVATOR2 Python/C++ port

Preparation and analysis use the [original automatic normalization method](docs/analysis-normalization.md):
size/MAP/GC correction, class-wise zero replacement, normalized-control pooling
and separate IN/OUT median centering. No calibration file or exposure estimates
are required.

[![Python port CI](https://github.com/danilotat/excavator2-py/actions/workflows/python-port.yml/badge.svg?branch=dev%2Fporting)](https://github.com/danilotat/excavator2-py/actions/workflows/python-port.yml)

This project reimplements **EXCAVATOR2**, a copy-number variant caller that uses
both in-target and off-target reads from whole-exome sequencing. Python owns the
workflow, reference and BAM access, normalization and reporting; a small C++ core
accelerates HSLM segmentation and FastCall numerical calculations.

`main` preserves the original scientific results, including observed legacy quirks.
On `dev/improvements`, scientific corrections are introduced independently.
FastCall corrections keep [fitting and reporting consistent](docs/fastcall-a19.md),
[optimize the truncated objective, and stabilize probabilities](docs/fastcall-fitting.md).
The original implementation is preserved
in `excavator2/`; the port lives in `src/excavator2/` and `cpp/`.

The three-stage pipeline on `main` was qualified against the original on the
supplied dataset and a differential test corpus. **Plotting is deferred.**
See the [compatibility report](docs/compatibility-qualification.md),
[performance results](docs/performance.md) and [roadmap](docs/python-cpp-porting-roadmap.md).

## CI performance comparison

The table is updated automatically from the latest passing benchmark on
`dev/porting`. Both versions run on the same Linux CI machine with matching
scientific outputs. This is a **small synthetic workload, not a real-WES or
full-pipeline speedup**. Target timing is one run; paired analysis uses three
measured runs with alternating execution order and one excluded warm-up.
Interpreter startup and artifact I/O are included; Docker setup and correctness
checks are excluded. Preparation timing is omitted because the original generates
plots that the port does not. Neither timed stage generates plots.

<!-- ci-benchmark:start -->

Latest passing [CI measurement](https://github.com/danilotat/excavator2-py/actions/runs/36269024286) · revision `a78d1dc` · 4,113 synthetic windows · one analysis test sample.

| Stage | Original | Python/C++ | Original / port | Measured runs |
| --- | ---: | ---: | ---: | ---: |
| Target generation | 26.860 s | 0.294 s | 91.32× | 1 |
| Paired analysis | 0.925 s | 0.389 s | 2.38× | 3 |

Analysis values are medians after one excluded warm-up. Ratios above 1 mean the port was faster.

<!-- ci-benchmark:end -->

Full timing and provenance data are retained in CI's `version-benchmark` artifact;
see [benchmark instructions](benchmarks/README.md) for local reproduction. Hosted
runner timings vary, so performance numbers are informational, not pass/fail thresholds.

## Installation

Requirements: **Python 3.11–3.13**, [uv](https://docs.astral.sh/uv/), and a C++17
compiler. Linux and macOS are tested. On macOS, install Apple Command Line Tools
with `xcode-select --install`; on Ubuntu/Debian, install `build-essential`.

```sh
git clone --branch dev/porting https://github.com/danilotat/excavator2-py.git
cd excavator2-py
uv sync --locked
uv run excavator2 --help
```

`uv sync` creates `.venv`, selects the repository's Python 3.12 version, installs
the locked dependencies and builds the C++ extension. The port does **not** need
R, Perl, Fortran or Docker at runtime. Docker is used only for original-version
comparisons below. This is a development package, not a published stable release.

To check the installation:

```sh
uv run pytest
```

See [development instructions](docs/development.md) for alternative installation,
build and testing commands.

## Inputs

Run commands from the repository root. Paths inside YAML files resolve relative
to the **working directory**, not the YAML file. Use the same assembly and
chromosome naming across all inputs.

Create `config.yaml` with your reference assets and capture BED:

```yaml
Reference:
  Assembly: hg38
  FASTA: data/reference.fa
  BigWig: data/mappability.bw
  Chromosomes: data/chromosome_coordinates.txt
  Centromeres: data/centromeres.txt
  Gaps: data/gaps.txt
Target:
  Name: MyCapture
  BED: data/capture.bed
  Window: 30000
```

The FASTA must be indexed. For example:

```sh
uv run python -c "import pysam; pysam.faidx('data/reference.fa')"
```

The chromosome-coordinate, centromere and gap tables use the legacy formats;
see [the example config](.test/config.yaml) and [reference annotation files](.test/ref/).
Reference FASTA and BigWig files are not bundled. The current compatibility policy
uses chromosomes 1–22 and X and rejects BED rows on unsupported contigs
(including Y/MT) instead of silently discarding them. Arbitrary contigs/assemblies
are not qualified. Gap-overlapping windows are excluded in full. Analysis rejects
windows touching or spanning a centromere; either chromosome arm may be absent.

Create `samples.yaml`, mapping unique sample names to coordinate-sorted,
indexed **BAM** files:

```yaml
Test1: data/test.bam
Control1: data/control.bam
```

Create `design.yaml`, mapping experimental labels to those same sample names:

```yaml
T1: Test1
C1: Control1
```

For paired analysis, each `T<number>` requires a matching `C<number>`. For pooling,
list the test samples as `T1`, `T2`, … and all controls as `C1`, `C2`, …; each test
uses the average of the normalized controls, following the original code.

## Run the port

Use a **new output directory** for each stage:

```sh
uv run excavator2 target --config config.yaml --output target/

uv run excavator2 prepare --samples samples.yaml --target target/ \
  --output prepared/ --threads 4

uv run excavator2 analyze --samples design.yaml --target target/ \
  --input prepared/ --output results/ --experiment paired \
  --parameters excavator2/parameters.yaml
```

Use `--experiment pooling` for pooled controls. Analysis currently uses one worker;
`--threads` in preparation bounds sample workers. Existing output directories and
`--force` are rejected. CRAM and plotting are not implemented/qualified.

The target directory contains a manifest and reference arrays; preparation saves
counts and normalization checkpoints. Under `results/Results/<sample>/`, analysis
writes:

- `HSLMResults_<sample>.txt`: window ratios and segmentation.
- `FastCallResults_<sample>.txt`: called regions and copy-number estimates.
- `EXCAVATORWindowCall_<sample>.vcf` and `EXCAVATORRegionCall_<sample>.vcf`.
- `checkpoints.npz`: numerical intermediates used for comparison.

The result manifest records inputs, parameters and execution settings. Random
FastCall ties require the original per-sample R RNG state via `--r-seed-states`;
see [the replay contract](docs/hslm-analysis.md#replaying-original-fastcall-random-ties).
Final probabilities use truncated class support, recorded in the manifest.
Probabilities of 1 inside a single allowed class express a hard model boundary,
not calibrated biological certainty. See [A19 and its remaining limits](docs/fastcall-a19.md).
Each command provides `--help`.

## Citation and license

Please cite D'Aurizio et al., *Enhanced copy number variants detection from
whole-exome sequencing data using EXCAVATOR2*, Nucleic Acids Research (2016),
44(20):e154, [doi:10.1093/nar/gkw695](https://doi.org/10.1093/nar/gkw695).
The original attribution and [CC BY-NC-SA 3.0 license](LICENSE) are retained.

## Run the original version and compare results

### Automated full-pipeline comparison

Start Docker with Linux/amd64 support, including emulation on Apple Silicon.
Pull the exact legacy image recorded by the project, then run the comparison:

```sh
legacy_image=$(uv run python -c "import json; print(json.load(open('tools/oracle/baseline.json'))['image'])")
docker pull --platform linux/amd64 "$legacy_image"

uv run python tools/oracle/end_to_end.py --output .oracle/comparison
```

This generates small reference files and BAMs, runs **all three stages of both
implementations**, and compares paired and pooling results. No human-reference
download is needed. It extracts the pinned original source and runs it in the
legacy container, using original native binaries only after verifying their
source identity. Both output trees, logs and `concordance.json` are retained.
Choose a new directory for every run.

### Run the original on the supplied WES example

Download the FASTA and BigWig from `reference_urls` in
[baseline.json](tools/oracle/baseline.json) to the paths in
[.test/config.yaml](.test/config.yaml), decompress the FASTA, and index it:

```sh
uv run python -c "import pysam, yaml; c=yaml.safe_load(open('.test/config.yaml')); pysam.faidx(c['Reference']['FASTA'])"

uv run python tools/oracle/oracle.py run --output .oracle/original --threads 1

uv run excavator2 target --config .test/config.yaml --output .oracle/new-target
uv run excavator2 prepare --samples .test/sample_sheet.yaml \
  --target .oracle/new-target --output .oracle/new-prepared --threads 4
uv run excavator2 analyze --samples .test/sample_file_list.yaml \
  --target .oracle/new-target --input .oracle/new-prepared \
  --output .oracle/new-results --experiment paired

uv run python tools/oracle/compare_analysis.py \
  .oracle/original/results .oracle/new-results --report .oracle/analysis-comparison.json
```

`oracle.py run` is deliberately restricted to the supplied, checksum-pinned paired
fixture; it is not a general runner for custom datasets. It invokes the original
`TargetPerla.pl`, `EXCAVATORDataPrepare.pl` and `EXCAVATORDataAnalysis.pl`, retaining
their outputs and exporting RData checkpoints. A normal clone is required because
the runner extracts the pinned source commit from Git history.

### Compare existing results, including your own data

Run the original and the port with identical BAMs, references, targets, parameters
and experimental design, keeping their output directories separate. Use `main`
when checking exact legacy equivalence; improvements may deliberately differ. Once both
runs have finished:

```sh
uv run python tools/oracle/compare_analysis.py \
  /path/to/original-results /path/to/rewritten-results --report comparison.json
```

Pass the directories **containing** `Results/`, not individual sample folders.
The comparator checks output-file inventory, HSLM coordinates/order/classes and
segment boundaries, FastCall tables and both VCFs. Calls must match exactly;
continuous HSLM values use the established tolerances. Only VCF `fileDate` is
ignored. A mismatch exits nonzero; a successful comparison writes a JSON report.
It does not compute reciprocal-overlap scores or chromosome-level summary charts.

For target/preparation checkpoint comparisons, use the automated full-pipeline
comparison above or the export/conversion tools described in the
[legacy-oracle guide](tools/oracle/README.md). Plots are excluded. Check original
logs and required output files as well: the original wrappers can report success
after an internal failure. The supplied shallow example has no non-normal calls;
the generated paired/pooling comparison also checks explicit losses and gains.
