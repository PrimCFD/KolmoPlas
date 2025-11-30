# `src/numerics/` — Finite-Volume Numerics Kernels (WIP)

The `numerics` directory contains a small, self-contained library of **finite-volume math kernels** used across physics modules. It focuses on **MAC-grid** (staggered) operators and heavy numerical loops, implemented in Fortran with thin C and C++ wrappers.

Numerics is deliberately:

* **Stateless** — no global simulation state.
* **CFD-agnostic** — no knowledge of specific equations/physics.
* **Interop-friendly** — C ABI + type-safe C++ API.

Higher-level orchestration (meshes, fields, I/O, MPI, etc.) lives outside this directory.

---

## 1. Directory layout

```text
numerics/
├─ CMakeLists.txt          # Build + optimisation rules for the numerics static library
├─ include/
│  └─ numerics/
│     └─ MacOps.hpp        # Public C++ API (namespace numerics::kernels)
├─ src/
│  └─ MacOps.cpp           # C++ implementations / wrappers around C/Fortran kernels
└─ kernels/
   ├─ kernels.f90          # Fortran implementations of MAC-grid kernels
   └─ kernels.h            # C ABI header (extern "C" prototypes)
```

This library is usually built as a static library named **`numerics`** and linked into higher-level components.

---

## 2. Responsibilities

Numerics **does**:

* Implement **low-level operators** on a MAC grid:

  * Divergence, gradients, and pressure-correction helpers.
  * Explicit/implicit diffusion schemes.
  * Nonlinear advection (e.g. KK3-style schemes).
  * Sub-grid-scale viscosity (e.g. Smagorinsky).
* Provide a stable **C ABI** via `kernels.h`.
* Wrap the C ABI in a simple, type-safe **C++ interface** under `numerics::kernels`.

Numerics **does not**:

* Own meshes, fields, or boundary conditions.
* Perform halo exchanges or MPI communication.
* Handle I/O, logging, or configuration.

Callers must prepare data (ghost cells, BCs, MPI halos, etc.) and then hand off raw pointers and sizes to the kernels.

---

## 3. Data layout & conventions

All kernels follow the same core conventions:

* **Structured Cartesian mesh** with a uniform ghost width `ng`.
* All sizes passed to kernels are **ghost-inclusive** totals, e.g.

  * `nxu_tot`, `nyu_tot`, `nzu_tot` for face-centred `u`.
  * `nxc_tot`, `nyc_tot`, `nzc_tot` for cell-centred scalars.
* Arrays are stored in a **3D, i-fastest** layout (row-major indexing).
* Velocity components `u`, `v`, `w` are **face-centred (MAC grid)**.
* Scalars (`p`, `rho`, `nu_eff`, `nu_t_c`, etc.) are **cell-centred**.

Callers are responsible for:

* Allocating arrays with the correct total sizes (including ghosts).
* Keeping ghost layers up-to-date before invoking kernels.
* Passing consistent `dx`, `dy`, `dz`, and time step `dt`.

---

## 4. C++ API overview (`numerics::kernels`)

The main user-facing header is:

```cpp
#include <numerics/MacOps.hpp>
```

Functions live in the `numerics::kernels` namespace and forward to C/Fortran kernels.

### 4.1 Typical kernel groups

At a high level, the API is organised into the following groups (each can be expanded with full signatures and examples):

1. **Sub-grid-scale viscosity**

   * Compute Smagorinsky eddy viscosity `nu_t` at cell centres from MAC velocities.

2. **Divergence and pressure gradient**

   * `∇·u` from MAC velocities to cell centres.
   * Face-centred pressure gradients from cell-centred pressure.

3. **Velocity correction / projection**

   * Velocity update with constant density.
   * Velocity update with variable density (cell-centred `rho_c`).

4. **Viscous diffusion**

   * Explicit forward-Euler step for diffusing MAC velocities.
   * Implicit backward-Euler sweeps (Jacobi and red–black Gauss–Seidel).
   * Residual evaluation for convergence monitoring.

5. **MAC advection**

   * Nonlinear advection terms for MAC velocities (e.g. KK3 scheme).

---

## 5. C ABI and Fortran layer

Internally, the heavy loops live in `kernels.f90` and are exposed via a C ABI declared in `kernels/kernels.h`.

* All kernels use **plain C types**: pointers to `double` and `int`-like sizes.
* Fortran routines are marked with `bind(C)` to keep symbol names predictable.
* `MacOps.cpp` implements the C++ wrappers, handling any small API niceties (e.g. mapping pointer output arguments to references).

The C ABI is available for integration with C or external solvers that do not want the C++ layer.

---

## 6. Build & optimisation (CMake)

`numerics/CMakeLists.txt` defines the `numerics` library and configures it for high performance:

* Compiles Fortran kernels with aggressive optimisation flags (e.g. `-O3` / `-Ofast`, vectorisation, loop unrolling), with per-compiler tweaks.
* Requires OpenMP and links against both the C++ and Fortran OpenMP runtimes.
* Optionally enables vectorisation reports (via a CMake option) to inspect compiler behaviour.

The library itself has **no public dependency** on higher-level simulation modules; it is intended to be reusable.

---

## 7. Typical usage pattern

From a physics or solver module, usage typically looks like:

1. Own velocity and scalar fields (with ghosts) in a higher-level data structure.
2. Apply boundary conditions and halo exchanges.
3. Extract raw pointers and ghost-inclusive extents.
4. Call an appropriate kernel, e.g. for divergence, diffusion, or advection.
5. Use the results to build RHS terms, perform projection, etc.

---

## 8. Adding a new kernel

When introducing a new numerical kernel to this directory:

1. **Implement** it in `kernels.f90` as a Fortran routine with `bind(C)` and a simple C-like signature.
2. **Declare** it in `kernels/kernels.h` inside the `extern "C"` block.
3. **Wrap** it in C++:

   * Add a declaration in `include/numerics/MacOps.hpp` under `namespace numerics::kernels`.
   * Implement the wrapper in `src/MacOps.cpp` calling the C symbol.
4. Keep the kernel side-effect-free except for its output arrays.
5. Follow the standard data layout and ghost-cell conventions.
