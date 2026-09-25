# Developing the Python-first port

This is an installable development package, version `0.1.0.dev0`, with a tested
[FastCall implementation](fastcall-reference.md). The CLI now generates targets,
prepares indexed BAMs and analyzes their counts. The all-new pipeline matches the
original calls on the supplied paired dataset; see [M5 usage and limits](target-porting.md).
The original implementation remains in `excavator2/`, and
its installation/workflow instructions remain in the root README.

## Prerequisites

- CPython 3.11–3.13; `.python-version` selects 3.12 for development.
- A C++17 compiler: Apple Command Line Tools on macOS or GCC/Clang and Python
  development headers on Linux.
- `uv` for the locked development environment (recommended), or Python/pip.
- Internet access for the first dependency install. Build isolation supplies
  CMake/Ninja when required; no global CMake installation is necessary.

Python dependencies: NumPy, SciPy, PyYAML, pysam, and pyBigWig. Plotting is an
optional extra. Our own native extension contains scalar FastCall and HSLM kernels. R, Perl, Fortran and samtools executables
are needed only for legacy-oracle regeneration. Target generation requires the
FASTA, BigWig and coordinate annotations referenced by the input YAML.

## Local setup

From the repository root:

```sh
uv sync --locked
uv run excavator2 --help
uv run excavator2 target --help
uv run python -c "from excavator2 import _core; print(_core.fastcall_posterior.__name__)"
uv run pytest
uv run ruff check .
uv run ruff format --check .
```

`uv sync` creates `.venv/` and installs the package in editable mode. Python edits
are visible immediately; rerun `uv sync --reinstall-package excavator2` after
changing C++ or CMake files. Use `uv sync --locked --extra plots` when needed.
The committed `uv.lock` fixes Python dependency resolution; build-system versions
are pinned separately in `pyproject.toml`. This is not yet a fully pinned compiler
or legacy runtime environment. Update locks deliberately with `uv lock`.

Alternatively, with a supported Python interpreter:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -e '.[plots]' pytest ruff
.venv/bin/python -m pytest
```

The pip route respects dependency bounds but does not consume `uv.lock`.

## Folder responsibilities

| Location | Purpose |
|---|---|
| `src/excavator2/cli.py` | Three-stage CLI and options |
| `src/excavator2/reference/` | Readable Python versions of accelerated kernels |
| `cpp/bindings/` | Minimal Python/C++ boundary |
| `cpp/include/excavator2/`, `cpp/tests/` | Kernel declarations and sanitizer smoke tests |
| `tests/` | Packaging, CLI, and differential numerical tests |
| `tests/differential/`, `tests/fixtures/` | Oracle comparisons and small goldens |
| `tools/oracle/` | Pinned legacy baseline and fixture generation |
| `benchmarks/` | Reproducible FastCall performance check |
| `excavator2/`, `.test/` | Preserved original implementation and example inputs |

`target_pipeline.py` owns target orchestration, `target.py` geometry, and
`target_features.py` reference extraction. Preparation, normalization, HSLM and
FastCall policy remain in their Python stage modules. Keep algorithm policy in Python. Custom C++ starts with
HSLM and FastCall hot kernels; further additions require evidence.

## Build and validation

```sh
uv build
```

This builds a source distribution and a wheel under `dist/`; the wheel includes
`excavator2._core`. There is no AVX requirement or host-specific ISA flag. Native
floating-point flags disable fast-math and contraction to preserve numerical compatibility. Wheels remain specific to their Python/platform ABI.

The separate `python-port.yml` workflow checks Python 3.11–3.13 on Linux and macOS,
including an installed-wheel test outside the checkout. This matrix is a target
for CI qualification, not a claim that every platform has already passed.
The original legacy workflow is unchanged. Differential tests cover native/Python
kernels and captured original geometry, features, preparation and analysis.
Live CI regenerates original geometry/features and preparation/analysis cases.
The target command now writes versioned artifacts consumed by prepare/analyze;
raw RData is not read directly. Local M5 supplied-data acceptance and 137 tests
passed in editable and installed-wheel environments. Broader M6 qualification
remains, including confirmation of the latest integration on remote CI.
See [the baseline report](legacy-baseline.md), [M3](hslm-analysis.md),
[M4](preparation.md) and [M5](target-porting.md).
See [the detailed roadmap](python-cpp-porting-roadmap.md).

Build follows the official [scikit-build-core guide](https://scikit-build-core.readthedocs.io/en/latest/guide/getting_started.html)
and [pybind11 CMake documentation](https://pybind11.readthedocs.io/en/stable/compiling.html).
Existing attribution and the root LICENSE are retained.
