# Developing the Python-first port

This is an installable development package, version `0.1.0.dev0`, with a tested
[FastCall implementation](fastcall-reference.md). The CLI does **not** yet prepare
samples or call CNVs. The original implementation remains in `excavator2/`, and
its installation/workflow instructions remain in the root README.

## Prerequisites

- CPython 3.11–3.13; `.python-version` selects 3.12 for development.
- A C++17 compiler: Apple Command Line Tools on macOS or GCC/Clang and Python
  development headers on Linux.
- `uv` for the locked development environment (recommended), or Python/pip.
- Internet access for the first dependency install. Build isolation supplies
  CMake/Ninja when required; no global CMake installation is necessary.

Python dependencies: NumPy, SciPy, PyYAML, pysam, and pyBigWig. Plotting is an
optional extra. Our own native extension contains three scalar FastCall batch kernels. R, Perl, Fortran, samtools executables, and full reference
assets are needed only for the legacy-oracle regeneration, not this package setup.

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

`fastcall.py` now owns the first numerical port. Add other Python stage modules
(`target.py`, `reads.py`, `normalization.py`, `design.py`, `hslm.py`) when implementing
those stages, rather than populating
empty placeholder APIs. Keep algorithm policy in Python. Custom C++ starts with
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
The original legacy workflow is unchanged. Scaffold tests check packaging and safe CLI failure. Separate FastCall differential
tests check both native and Python backends against original-R fixtures; end-to-end scientific
equivalence remains unproven.

All three stage commands accept their planned options but exit nonzero without
writing output until implemented. Continue using the original scripts for real
analysis. A repeatable original-pipeline baseline has now been captured; see
[the baseline report](legacy-baseline.md). M2 is complete: both FastCall backends match the saved fixtures. Python retains
initialization, iteration/stopping, and assignment. M3 has not started.
See [the detailed roadmap](python-cpp-porting-roadmap.md).

Build follows the official [scikit-build-core guide](https://scikit-build-core.readthedocs.io/en/latest/guide/getting_started.html)
and [pybind11 CMake documentation](https://pybind11.readthedocs.io/en/stable/compiling.html).
Existing attribution and the root LICENSE are retained.
