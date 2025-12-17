<p align="center">
  <img src="assets/logo/KolmoPlasLogo512.png" alt="KolmoPlas logo" width="160">
</p>

# KolmoPlas

> **Status:** _Hydrodynamics_

[![Docs](https://img.shields.io/badge/docs-online-blue)](https://primcfd.github.io/KolmoPlas/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Linux CI](https://github.com/PrimCFD/KolmoPlas/actions/workflows/linux.yml/badge.svg?branch=main)](https://github.com/PrimCFD/KolmoPlas/actions/workflows/linux.yml)
[![Docs](https://github.com/PrimCFD/KolmoPlas/actions/workflows/docs.yml/badge.svg?branch=main)](https://github.com/PrimCFD/KolmoPlas/actions/workflows/docs.yml)
[![Style](https://github.com/PrimCFD/KolmoPlas/actions/workflows/style.yml/badge.svg?branch=main)](https://github.com/PrimCFD/KolmoPlas/actions/workflows/style.yml)


---

## 1&nbsp;· Overview

The project is an **open, modular finite‑volume 3D LES solver for thermal–plasma flows** that lets researchers drop in new sub‑grid‑scale (SGS) models with minimal friction while still scaling to modern CPU & GPU clusters. The project is intended to research new SGS models for thermal plasma jets simulations accounting for steep density and temperature gradients.

---

## 2&nbsp;· Repository Layout (high‑level)

```text
KolmoPlas/
├─ .github/workflows/        # CI definitions
├─ assets/                   # Docs and visual assets
├─ cmake/                    # CMake helper modules
├─ docs/                     # Sphinx/Doxygen furo doc
├─ examples/                 # Tiny example run‑ready cases
├─ extern/                   # Third‑party sources (tar balls mirrored, no edits)
├─ scripts/                  # Dev‑ops & helper scripts
├─ src/                      # Solver source code
│   ├─ apps/                 # Executable entry points 
│   ├─ core/                 # C++ runtime static lib (orchestration, memory management)
│   ├─ gui/                  # Qt/VTK front‑end
│   ├─ ipc/                  # Inter-process communication GUI/Solver
│   ├─ numerics/             # Common math kernels static lib
│   └─ physics/              # Hot‑swappable physics modules dynamic lib
├─ tests/                    # Unit, regression & perf tests
├─ .clang-format             # C++ style
├─ .cmake-format.yml         # CMake style
├─ .fprettify.yml            # Fortran style
├─ .gitignore
├─ CMakeLists.txt            # Root super‑build
├─ kolmoplas                 # CLI mini-tool
├─ LICENSE
└─ README.md                 # This file
```

---

## 3&nbsp;· Directory Details

| Path | Role | Notes |
|------|------|-------|
| **src/core** | Owns mesh, fields, time loop, I/O; _no physics_. | Written in modern C++ (C++20). |
| **src/physics** | Physics shared library plugins | Each physics → one shared lib (`libphysics_*.so`). |
| **src/numerics** | Optimised Fortran math; reused by plug‑ins. | Vectorised / GPU‑offloaded via OpenMP/CUDA Fortran. |
| **extern** | Vendored: CGNS, PETSc, etc. | Pulled via `FetchContent`/`ExternalProject`; do **not** modify in‑tree. |
| **tests** | Unit, regression CGNS, etc | CI or manually ran. |
| **examples** | Minimal decks: Taylor–Green, etc. | Small test run and CFD/MHD benchmarks. |

---

## 4&nbsp;· Quick Start

This project ships convenience scripts in `scripts/` for reliable, repeatable developer workflows (builds, docs, MPI, offline vendor cache, cleaning, formatting). See the cheatsheet in that folder for details.

The scripts are accessible through the thin `kolmoplas` CLI tool.

### 4.1 Prerequisites
- **CMake ≥ 3.24**, a C/C++/Fortran toolchain; **Ninja** is recommended.
- **Parrallel stack** (OpenMP) technically overridable but no gains.
- Optional/when needed:
  - **MPI stack** (Open MPI) for MPI builds and tests.
  - **Doxygen** and **Sphinx** (`sphinx-build`) for docs.
  - Formatters: `clang-format`, `fprettify`, `cmake-format`. (if using a Python virtual environment, run `source venv/bin/activate`)


> Tip: CI uses recent gcc on Linux, might be tweaking required for clang; matching that locally avoids surprises (see §5 CI). System gcc + OpenMP typically is the HPC setup.

### 4.2 `kolmoplas` CLI overview

The `kolmoplas` helper is a **thin Python CLI** that simply routes to the shell/Python scripts in `scripts/`. Its goals are:

* ergonomic, memorable commands (`./kolmoplas build`, `./kolmoplas test reg`, …)
* stable entry points for docs/CI, even if internal scripts move or grow
* a single place to discover the developer workflow via `-h/--help`.

Basic usage:

```bash
# Top-level help
./kolmoplas -h

# Help for a specific subcommand or group
./kolmoplas build -h
./kolmoplas test -h
./kolmoplas post -h
```

Top-level commands currently map as follows:

* `build` → `scripts/build.sh`
  Configure & build (`BUILD_DIR`, `CMAKE_BUILD_TYPE`, etc. are passed through).

* `clean` → `scripts/clean_build.sh`
  Clean CMake/Ninja build trees (e.g. `./kolmoplas clean --all -y`).

* `docs` → `scripts/build_docs.sh`
  Build and optionally serve Sphinx/Doxygen docs (`--serve`, `--open`, …).

* `test` → family of test runners:

  * `./kolmoplas test unit` → `scripts/run_unit_tests.sh`
  * `./kolmoplas test perf` → `scripts/run_perf_tests.sh`
  * `./kolmoplas test reg` → `scripts/run_regression_tests.sh`

* `scaling` → scaling studies:

  * `./kolmoplas scaling run`
    wraps `scripts/run_scaling.sh` (sweeps grids × MPI ranks × threads and writes `build-scaling/results_scaling.csv`).
  * `./kolmoplas scaling post`
    wraps `scripts/postprocess_scaling.py` to clean the CSV and generate speedup/efficiency plots.

* `clean-extern` → `scripts/clean_extern.sh`
  Wipe `extern/` vendor caches (keeps `extern/README.md`).

* `fmt` → `scripts/format_all.sh`
  Run C/C++/Fortran/CMake formatters (honours `STRICT=1`).

* `prefetch` → `scripts/prefetch_third_party.sh`
  Populate `extern/` with reproducible third-party archives for offline builds.

* `post` → post-processing utilities:

  * `./kolmoplas post tgv-spectrum`
    → `scripts/tgv_spectrum.py` (isotropic TGV energy spectrum from CGNS/HDF5, saves plots and data near the input file).
  * `./kolmoplas post ghia`
    → `scripts/ldc_ghia.py` (lid-driven cavity centerline profiles vs Ghia et al., errors + plots).

* `examples` → regression-backed solver examples:

  * `./kolmoplas examples run <case>`
    runs a named regression case (e.g. `tgv128`, `channel128`) by forwarding a CTest name regex into `scripts/run_regression_tests.sh`.

Debugging tip: set `KOLMOPLAS_DEBUG=1` to have the CLI echo the underlying script command before running it.

```bash
KOLMOPLAS_DEBUG=1 ./kolmoplas test reg
```

### 4.3 Typical developer workflow with `kolmoplas`

A common day-to-day workflow using the CLI looks like:

1. **One-time (or infrequent) vendor cache**

   ```bash
   # Mirror third-party archives into extern/
   ./kolmoplas prefetch
   ```

   Later you can build fully offline with:

   ```bash
   ./kolmoplas build OFFLINE=1
   ```

2. **Build the solver**

   ```bash
   # Default Release build
   ./kolmoplas build

   # Debug build in its own tree
   ./kolmoplas build BUILD_DIR=build-debug CMAKE_BUILD_TYPE=Debug
   ```

3. **Run tests**

   ```bash
   # Fast unit tests
   ./kolmoplas test unit

   # Performance/sanity runs
   ./kolmoplas test perf

   # MPI regression suite (laptop- and cluster-friendly)
   ./kolmoplas test reg
   ```

4. **Run an example case**

   ```bash
   # Backed by regression tests; case name comes from CTest
   ./kolmoplas examples run tgv128

   # Or a cavity/channel case, when present
   ./kolmoplas examples run cavity256
   ```

5. **Post-process results**

   ```bash
   # TGV energy spectrum from the latest TGV regression output
   ./kolmoplas post tgv-spectrum --data-dir build-regression/tests/regression/tgv128/out

   # Lid-driven cavity centerline comparison vs Ghia et al.
   ./kolmoplas post ghia --data-dir build-regression/tests/regression/cavity256/out
   ```

6. **Scaling studies**

   ```bash
   # Sweep grids × ranks × threads and log to build-scaling/
   ./kolmoplas scaling run

   # Produce cleaned CSV + strong-scaling & efficiency plots
   ./kolmoplas scaling post
   ```

7. **Docs + housekeeping**

   ```bash
   # Build or preview docs
   ./kolmoplas docs
   ./kolmoplas docs --serve --open

   # Clean build trees and (optionally) vendor caches
   ./kolmoplas clean --all -y
   ./kolmoplas clean-extern -y

   # Run formatters
   ./kolmoplas fmt
   ./kolmoplas fmt STRICT=1
   ```

---

## 5&nbsp;· Continuous Integration
[![Linux CI](https://github.com/PrimCFD/KolmoPlas/actions/workflows/linux.yml/badge.svg?branch=main)](https://github.com/PrimCFD/KolmoPlas/actions/workflows/linux.yml) [![Docs](https://github.com/PrimCFD/KolmoPlas/actions/workflows/docs.yml/badge.svg?branch=main)](https://github.com/PrimCFD/KolmoPlas/actions/workflows/docs.yml) [![Style](https://github.com/PrimCFD/KolmoPlas/actions/workflows/style.yml/badge.svg?branch=main)](https://github.com/PrimCFD/KolmoPlas/actions/workflows/style.yml)

All CI is always built from clean slate (containerized) to check the whole pipeline on Hardware with missing tools (run on GitHub manually)

* **linux.yml** – GCC 13 / Clang 18 matrix; runs unit and performance.
* **style.yml** – clang‑format, fprettify, cmake‑lint. (only CI automatically ran on push)
* **docs.yml** – builds Sphinx docs, pushes to `gh-pages`.

---

## 6&nbsp;· License

Distributed under the **MIT License** as outlined in `LICENSE`.