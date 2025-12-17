# `scripts/` Developer Workflow Cheatsheet

## TL;DR table

| Script | Purpose | Typical use |
|---|---|---|
| `build.sh` | Configure & build the project via CMake (env‑driven). | `./scripts/build.sh` |
| `build_docs.sh` | Build clean Sphinx HTML (+ Doxygen XML) and optionally serve locally. | `./scripts/build_docs.sh --serve --open` |
| `clean_build.sh` | Remove CMake/Ninja build trees (`build-*`, keeps docs). | `./scripts/clean_build.sh --all -y` |
| `clean_extern.sh` | Wipe vendored third‑party sources under `extern/` (keeps `extern/README.md`). | `./scripts/clean_extern.sh -y` |
| `format_all.sh` | Run code formatters for C/C++, Fortran, and CMake (skips vendor/build). | `./scripts/format_all.sh` |
| `mpi_env.sh` | Sourceable helpers to set MPI launcher/env; provides `mpi_exec`. | `source scripts/mpi_env.sh auto` |
| `prefetch_third_party.sh` | Pre‑download third‑party sources and create reproducible archives in `extern/`. | `./scripts/prefetch_third_party.sh` |
| `run_unit_tests.sh` | Build (optional) and run **unit** tests via CTest; write JUnit XML. | `./scripts/run_unit_tests.sh` |
| `run_perf_tests.sh` | Build (optional) and run **perf** tests via CTest; write JUnit XML. | `./scripts/run_perf_tests.sh` |
| `run_regression_tests.sh` | Configure with MPI and run **MPI** integration tests, laptop/cluster friendly. | `./scripts/run_regression_tests.sh` |
| `run_scaling.sh` | Sweep grids × MPI ranks × OpenMP threads and write strong-scaling CSVs/logs. | `./scripts/run_scaling.sh` |
| `tgv_spectrum.py` | Compute isotropic TGV energy spectrum E(k) from CGNS/HDF5 outputs. | `python scripts/tgv_spectrum.py --file result.cgns` |
| `ldc_ghia.py` | Lid-driven cavity centerline regression vs Ghia et al. (profiles, errors, plots). | `python scripts/ldc_ghia.py --file cavity.cgns` |
| `postprocess_scaling.py` | Clean and plot strong-scaling CSVs (speedup/efficiency reports). | `python scripts/postprocess_scaling.py` |

> Unless stated otherwise, run scripts from the repo root. Many respect `BUILD_DIR` if you want custom build folders.

---

## Common prerequisites

- **CMake ≥ 3.24** and a C/C++/Fortran toolchain. Ninja recommended.
- **Python 3** (only for `build_docs.sh --serve`).
- Optional formatters: `clang-format`, `fprettify`, `cmake-format` (`STRICT=1` makes formatting failures fatal).
- For MPI workflows: an MPI stack (Open MPI/MPICH) and, on clusters, a launcher such as `srun`.
- For docs: `doxygen` and `sphinx-build` (the script will report if missing).

---

## `build.sh` — one‑stop CMake build

**What it does**  
- Picks a generator (prefers Ninja), configures CMake in `BUILD_DIR` and builds with parallel jobs.  
- Respects and forwards FetchContent & third‑party hints; can work **offline** if archives are present in `extern/`.  
- Filters noisy third‑party subbuild logs while still surfacing real errors.

**Key env vars (all optional)**  
- `BUILD_DIR` — build directory (default: `build`).  
- `CMAKE_BUILD_TYPE` — `Debug|Release|RelWithDebInfo|MinSizeRel` (default: `Release`).  
- `BUILD_TESTS` — `ON|OFF` (default: `ON`).  
- `MPIEXEC_PREFLAGS` — forwarded to CMake/CTest (e.g. `--oversubscribe`).  
- `CMAKE_GENERATOR` — override generator (e.g. `Ninja`).  
- `CMAKE_TOOLCHAIN_FILE` — forwarded as‑is.  
- `EXTRA_CMAKE_ARGS` — extra space‑separated CMake args.  
- `OFFLINE=1` — force disconnected `FetchContent` if you’ve already mirrored deps in `extern/`.
- `NPROCS` — override auto‑detected parallel job count.

**Examples**  
```bash
# Default Release build with Ninja
./scripts/build.sh

# Debug build into a custom dir
BUILD_DIR=build-debug CMAKE_BUILD_TYPE=Debug ./scripts/build.sh

# Pass OpenMPI oversubscribe
MPIEXEC_PREFLAGS=--oversubscribe ./scripts/build.sh

# Extra CMake toggles
EXTRA_CMAKE_ARGS="-DBUILD_EXAMPLES=ON -DBUILD_GUI=OFF" ./scripts/build.sh
```

---

## `build_docs.sh` — docs builder/previewer

**What it does**  
- Builds only the documentation: Doxygen XML → Sphinx HTML (no full CMake build).  
- Uses fixed paths under `build-docs/` to match `docs/conf.py`.  
- Can start a local HTTP server for quick preview.

**Flags**  
- `--clean` — remove previous docs build.  
- `--force` — clean + rebuild (CI‑friendly).  
- `--no-doxygen` — skip Doxygen step, reuse existing XML.  
- `--serve [--port N]` — serve HTML locally (default port 8000).  
- `--open` — try to open a browser tab after `--serve`.  
- `--quiet` — less verbose output.

**Examples**  
```bash
# Build docs once
./scripts/build_docs.sh

# Force clean + rebuild (useful in CI)
./scripts/build_docs.sh --force

# Rebuild and preview locally
./scripts/build_docs.sh --serve --open --port 8001
```

---

## `clean_build.sh` — clean build trees

**What it does**  
- Deletes CMake/Ninja artifacts in selected build directories.  
- Skips `build-docs/`

**Flags**  
- `-a, --all` — clean all `build-*` directories found.  
- `-k, --keep-deps` — preserve `_deps/` (downloaded third‑party sources) inside each build dir.  
- `-n, --dry-run` — show what would be removed.  
- `-y, --yes` — do not prompt.  
- `-v, --verbose` — more logging.  
- `-h, --help` — usage.

**Examples**  
```bash
# Clean a specific dir
./scripts/clean_build.sh build-debug

# Preview what would be deleted
./scripts/clean_build.sh --all -n

# Clean multiple builds but keep vendor downloads
./scripts/clean_build.sh --keep-deps build-regression build-perf -y
```

---

## `clean_extern.sh` — wipe vendored sources

**What it does**  
- Removes everything under `extern/` **except** `extern/README.md`.  
- Supports dry‑run and prompt‑less modes.

**Flags**  
- `-y, --yes` — do not prompt.  
- `-n, --dry-run` — list what would be removed.

**Examples**  
```bash
# See what would be deleted
./scripts/clean_extern.sh -n

# Nuke vendored archives/sources (keeps README)
./scripts/clean_extern.sh -y
```

---

## `format_all.sh` — code formatting

**What it does**  
- Runs the available formatters over project sources, skipping vendor/build trees:  
  - C/C++ via `clang-format`  
  - Fortran via `fprettify`  
  - CMake via `cmake-format`  
- Ignores missing tools unless `STRICT=1`.

**Env**  
- `STRICT=1` — fail the script if any formatter is missing or errors.

**Examples**  
```bash
# Best‑effort formatting
./scripts/format_all.sh

# CI‑style strict mode
STRICT=1 ./scripts/format_all.sh
```

---

## `mpi_env.sh` — sourceable MPI helpers

**What it does**  
- Detects MPI vendor/launcher and exports CMake/CTest‑friendly variables:  
  `MPIEXEC_EXECUTABLE`, `MPIEXEC_NUMPROC_FLAG`, `MPIEXEC_PREFLAGS`, `MPIEXEC_POSTFLAGS`.  
- Provides a convenience function `mpi_exec <np> <cmd …>` that calls the right launcher (`srun`, `mpirun`, `mpiexec`).  
- Modes:  
  - `auto` — detect cluster vs laptop, choose safe defaults.  
  - `emulate` — force laptop‑friendly settings (loopback/TCP/oversubscribe for Open MPI).  
  - `cluster` — disable emulation and favor scheduler/HCAs.

**Usage**  
```bash
# Source it (MUST be sourced)
source scripts/mpi_env.sh auto

# Then launch binaries consistently
mpi_exec 4 ./build/bin/your_mpi_program
```

---

## `prefetch_third_party.sh` — offline vendor cache

**What it does**  
- Configures CMake in a temporary dir with `-DPREFETCH_THIRD_PARTY=ON`, letting the project's CMake download declared third‑party sources.  
- Packages each fetched source as a reproducible archive into `extern/`, and writes `MANIFEST.prefetch` and `SHA256SUMS`.  
- Speeds up CI and enables fully offline builds when paired with `OFFLINE=1` in `build.sh`.

**Example**  
```bash
./scripts/prefetch_third_party.sh
# Archives now in extern/, plus MANIFEST.prefetch and SHA256SUMS
```

---

## Test runners

### `push_ci.sh`

**What it does**
- Pushes changes to GitHub  
- Configures tmux terminal with the actions user and actions-runner/run.sh  
- Waits on and runs CI config /.workflow/linux.yml  

**Example**  
```bash
./scripts/push_ci.sh
```

### `run_unit_tests.sh`
- **What it does:** Optionally builds (via `build.sh`) then runs `ctest -L unit`. Creates JUnit XML under `$BUILD_DIR/test-reports/unit/` (fallback created if CTest/Catch2 doesn’t emit one).  
- **Env:**  
  - `BUILD_DIR` (default: `build-unit`)  
  - `SKIP_BUILD` (`0|1`, default `0`)  
  - `CTEST_PARALLEL_LEVEL` (defaults to CPU count)  
  - `CTEST_TIMEOUT` (default: `900`)  
  - `REPORT_DIR` (default: `$BUILD_DIR/test-reports/unit`)  
- **Example:**  
  ```bash
  ./scripts/run_unit_tests.sh
  SKIP_BUILD=1 ./scripts/run_unit_tests.sh
  ```

### `run_perf_tests.sh`
- **What it does:** Same shape as unit, but runs `ctest -L perf` and writes to `$BUILD_DIR/test-reports/perf/`.  
- **Env:** `BUILD_DIR` (default `build-perf`), `SKIP_BUILD`, `CTEST_PARALLEL_LEVEL`, `CTEST_TIMEOUT`, `REPORT_DIR`.  
- **Example:**  
  ```bash
  CMAKE_BUILD_TYPE=Release ./scripts/run_perf_tests.sh
  ```

### `run_regression_tests.sh`
- **What it does:** Sources `mpi_env.sh`, configures with MPI and runs cases with the solver. Emits JUnit to `$BUILD_DIR/test-reports/regression/` (with a safe fallback if no launcher is found).  
- **Examples:**  
  ```bash
  ./scripts/run_regression_tests.sh
  BUILD_DIR=build-regression ./scripts/run_regression_tests.sh
  ```

## Scaling studies

### `run_scaling.sh`

**What it does**  
- Sweeps a set of cubic grid sizes `N` (i.e. `N^3`) across combinations of MPI ranks and OpenMP threads.  
- Runs the Poisson benchmark (`bench_poisson_fluids`) under MPI, collecting:  
  - mean/StdDev of per-step timings (µs) when available,  
  - wall-clock time,  
  - PETSc KSP iteration counts and convergence reasons.  
- Writes:  
  - `build-scaling/build_info.csv` with build/git/launcher info,  
  - `build-scaling/results_scaling.csv` with one row per (grid, ranks, threads, repetition),  
  - logs for each run under `build-scaling/logs/`.

**Key env**  
- `SCALING_GRIDS` — space-separated grid sizes `N` (default: `64 128 256`).  
- `SCALING_RANKS` — MPI ranks to test (default: `1 2 4`).  
- `SCALING_THREADS` — `OMP_NUM_THREADS` values to test (default: `1 2`).  
- `SCALING_REPS` — repetitions per configuration (default: `3`).  
- `SCALING_BUILD_FIRST` — `1|0`, whether to rebuild perf tests first (default: `1`).  
- `SCALING_SET_DA_PROCS` — `1|0`, auto-choose near-cubic PETSc DMDA decomposition (default: `1`).  
- `BUILD_DIR` — scaling build directory (default: `build-scaling`).  
- `PETSC_OPTIONS_EXTRA` — extra PETSc options appended to the solver run.

**Example**  
```bash
# Default sweep (64^3, 128^3, 256^3; ranks 1/2/4; threads 1/2)
./scripts/run_scaling.sh

# Custom grids/ranks/threads with extra PETSc options
SCALING_GRIDS="128 256" \
SCALING_RANKS="1 4 8" \
SCALING_THREADS="1 2" \
PETSC_OPTIONS_EXTRA="-pc_type mg -ksp_type pipecg" \
  ./scripts/run_scaling.sh
```

### `postprocess_scaling.py`

**What it does**

* Reads `build-scaling/results_scaling.csv` and:

  * normalises timings (prefers `mean_us`, falls back to `wall_s`),
  * annotates rows with a boolean `nice_factorization` flag when the DMDA layout is “nice” for the given `N` and `ranks`,
  * aggregates to a cleaned table with one row per `(grid_n, ranks)` including speedup and parallel efficiency vs the minimum-rank run.
* Writes:

  * `build-scaling/scaling_reports/results_scaling_cleaned.csv`,
  * speedup plots `speedup_N{N}.png` and efficiency plots `efficiency_N{N}.png` for each grid size.

**CLI**

```bash
python scripts/postprocess_scaling.py \
  --csv build-scaling/results_scaling.csv \
  --out-dir build-scaling/scaling_reports
```

**Arguments**

* `--csv` — input CSV path (default: `build-scaling/results_scaling.csv`).
* `--out-dir` — output directory for cleaned CSV and plots (default: `build-scaling/scaling_reports`).

## Post-processing utilities

### `tgv_spectrum.py`

**What it does**  
- Reads a Taylor–Green vortex (or similar) result file produced by the CGNS/HDF5 writers (`u_cell`, `v_cell`, `w_cell` cell-centred velocities).  
- Supports both backends:
  - CGNS: picks the last `FlowSolutionAtStepXXXX_Cell` by default.  
  - HDF5: picks the last `Step_xxxxxx` group (or an explicit index).  
- Computes an isotropic 1D energy spectrum `E(k)` by:
  1. subtracting mean velocity,  
  2. performing 3D FFTs,  
  3. forming kinetic energy density in k-space,  
  4. binning into spherical shells (log-spaced by default).  
- Produces spectrum plots and optionally a compensated spectrum, saving figures and data next to the result file.

**Typical usage**  
```bash
# From a specific CGNS file
python scripts/tgv_spectrum.py --file build-regression/tests/regression/tgv128/out/flow.cgns

# Let the script auto-discover the newest file under a directory
python scripts/tgv_spectrum.py --data-dir build-regression/tests/regression/tgv128/out
```

**Common options**

* `--file` — explicit `*.cgns` / `*.h5` / `*.hdf5` result file.
* `--data-dir` — directory to search for the latest result when `--file` is omitted.
* `--solution-name` — choose a specific `FlowSolution_t` in CGNS.
* `--step-index` — choose a specific `Step_xxxxxx` in HDF5 (default: last step).
* `--no-plot` — skip showing plots (still writes files).

---

### `ldc_ghia.py`

**What it does**

* Reads a 2D lid-driven cavity result (CGNS or HDF5):

  * CGNS: `u_cell`, `v_cell` from the last `*_Cell` `FlowSolution_t`,
  * HDF5: `u_cell`, `v_cell` from the last `Step_xxxxxx` group.
* Extracts:

  * `u(y)` along the vertical centerline `x = 0.5`,
  * `v(x)` along the horizontal centerline `y = 0.5`.
* Loads Ghia et al. benchmark data from ASCII tables (`ghiau.txt`, `ghiav.txt`).
* Interpolates your profiles to the Ghia sample points, prints L2 / L∞ errors, and overlays plots.
* Saves profiles, reference data, and interpolated values to text files alongside the input result, plus PNG/SVG/PDF figures.

**Defaults**

* Results directory:
  `build-regression/tests/regression/cavity256/out`.
* Ghia data directory: `scripts/ghia/` (`ghiau.txt`, `ghiav.txt`).

**Typical usage**

```bash
# Use the latest cavity result under the default regression output dir
python scripts/ldc_ghia.py

# Explicit file and Reynolds number
python scripts/ldc_ghia.py --file build-regression/tests/regression/cavity256/out/cavity.cgns --Re 1000

# Explicit Ghia data and output prefix
python scripts/ldc_ghia.py \
  --file cavity.h5 \
  --ghia-u scripts/ghia/ghiau.txt \
  --ghia-v scripts/ghia/ghiav.txt \
  --out-prefix build-regression/tests/regression/cavity256/out/cavity_ghia_Re1000
```

**Key options**

* `--file` / `--data-dir` — choose the result file or search directory.
* `--Re` — Reynolds number to pick the appropriate column from Ghia tables (default: `1000`).
* `--ghia-u`, `--ghia-v` — override paths to the Ghia reference tables.
* `--out-prefix` — prefix for saved data/plots (default: derived from the input filename).
* `--formats` — comma-separated figure formats (e.g. `png,svg,pdf`).
* `--no-plot` — don’t open GUI windows; just write files.

---

## Tips

- Use `BUILD_DIR=…` consistently to keep separate trees for `unit`, `perf`, and `regression`.  
- `clean_build.sh --keep-deps` is handy when you don’t want to re‑download third‑party code.  
- Pair `prefetch_third_party.sh` with `OFFLINE=1 ./scripts/build.sh` for air‑gapped or flaky‑network environments.
