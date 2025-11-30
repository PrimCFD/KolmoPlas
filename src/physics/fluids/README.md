# `physics/fluids/` — Incompressible Fluids Physics Plugin

The **fluids** directory implements an incompressible Navier–Stokes solver as a **physics plugin** for the core framework. It glues together:

* Mesh and field ownership from `core::master`.
* Low-level MAC-grid numerics from the `numerics` library.
* PETSc-based pressure solves for the projection step.

All logic is encapsulated behind the plugin/registry interface, so applications can enable the fluids solver by selecting the `fluids` program in their run configuration.

---

## 1. Directory layout

```text
physics/fluids/
├─ CMakeLists.txt           # Build rules for the fluids plugin (physics_fluids)
├─ include/
│  └─ fluids/
│     ├─ ApplyBCs.hpp
│     ├─ InitTG.hpp
│     ├─ InitChannel.hpp
│     ├─ SGS.hpp
│     ├─ MomentumPredictor.hpp
│     ├─ PressurePoisson.hpp
│     ├─ VelocityCorrector.hpp
│     ├─ ProjectionLoop.hpp
│     └─ Program.hpp        # Top-level program wiring
└─ src/
   ├─ Register.cpp          # Plugin registration entry point
   ├─ ApplyBCs.cpp          # Boundary-condition application for u,v,w,p
   ├─ InitTG.cpp            # Taylor–Green vortex initial condition
   ├─ InitChannel.cpp       # Plane channel-flow initial condition
   ├─ SGS.cpp               # Smagorinsky sub-grid-scale model
   ├─ MomentumPredictor.cpp # Explicit/implicit momentum predictor
   ├─ PressurePoisson.cpp   # PETSc-based pressure Poisson solve
   ├─ VelocityCorrector.cpp # Pressure-gradient velocity correction
   ├─ ProjectionLoop.cpp    # Orchestrates one projection time step
   └─ Program.cpp           # Assembles the full fluids program
```

The CMake target is a shared library named **`physics_fluids`**, linked against `core`, `numerics`, and `PETSC::petsc`.

---

## 2. Responsibilities & scope

The fluids plugin is responsible for:

* Advancing **incompressible flow** on a structured MAC grid.
* Providing reusable **actions** (init, SGS, predictor, pressure solve, corrector, projection loop).
* Handling boundary conditions at the physics level (in terms of u, v, w, p).
* Owning all PETSc configuration for the pressure Poisson solve.

It deliberately does **not**:

* Implement low-level difference operators (those live in `numerics/`).
* Own meshes or fields (they are obtained from `core::master::FieldCatalog` and `MeshTileView`).
* Perform I/O or high-level orchestration beyond a single timestep.

---

## 3. Plugin wiring (`Register.cpp`, `Program.cpp`)

### 3.1 Registration entry point

`src/Register.cpp` exposes the plugin to the global registry via an `extern "C"` hook, e.g.:

* Registers the **program** name: `"fluids"`.
* Registers individual **actions** under names like:

  * `"fluids.init_tg"`
  * `"fluids.init_channel"`
  * `"fluids.sgs"`
  * `"fluids.momentum_predictor"`
  * `"fluids.pressure_poisson"`
  * `"fluids.velocity_corrector"`
  * `"fluids.projection_loop"`

These names are what configuration files refer to when wiring custom action pipelines.

### 3.2 Fluids program

`Program.hpp` / `Program.cpp` define **`FluidsProgram`**, a `plugin::Program` that builds the standard fluids time integrator from a key–value configuration:

* Parses global parameters like grid spacing, fluid density, viscosity, SGS constants, etc.
* Parses boundary-condition specifications of the form `bc.<field>.<face> = type:value`.
* Instantiates sub-actions:

  * Initial condition (`InitTG` or `InitChannel`, depending on config).
  * SGS model (`SGS`).
  * Momentum predictor (`Predictor`).
  * Pressure Poisson solve (`Poisson`).
  * Velocity corrector (`Corrector`).
  * Projection loop (`ProjectionLoop`) that stitches the above into one timestep.

When this README is expanded, this section is a good place for a **config reference table** describing every supported key and its default.

---

## 4. Boundary conditions (`ApplyBCs`)

`ApplyBCs` translates high-level BC specs into concrete ghost-cell fills on the MAC grid:

* Consumes a `Params::bcs` table built by `Program` from keys like `bc.u.north = dirichlet:0.0`.
* Supports types:

  * `dirichlet`
  * `neumann`
  * `extrap` (extrapolate)
  * `mirror` / `symmetry` (velocity reflection)
  * `periodic` (usually handled by halo exchange; treated as a no-op here)
* Uses `mesh::apply_scalar_bc` for scalars and `mesh::apply_mirror_vector` for velocity mirrors.
* Respects domain periodicity from the mesh (per-axis periodic flags), so that periodic faces are not overwritten by BC tables.

Phase: **PreExchange** — it runs before halo exchanges so that BC-filled ghosts participate in MPI communication.

---

## 5. Initial conditions (`InitTG`, `InitChannel`)

### 5.1 `InitTG` — Taylor–Green vortex

Provides a canonical incompressible test case:

* Parameters (all optional, with sensible defaults):

  * `Lx`, `Ly`, `Lz` — domain lengths.
  * `U0` — velocity scale.
* Assumes u, v, w are MAC-staggered (`u` on x-faces, etc.) and uses mesh metadata as the single source of truth for:

  * Ghost width (`ng`).
  * Local and global interior cell counts.
  * Global offsets of the current tile.
* Fills u, v, w with the TG analytic field and leaves pressure zero by default.

Phase: **PreExchange** — fields are filled before any halo exchange.

### 5.2 `InitChannel` — plane channel flow

Sets up a plane channel-flow with optional perturbation:

* Parameters:

  * `Ly` — wall-normal extent.
  * `Ubulk` — target bulk velocity.
  * `pert` — amplitude of a small perturbation field.
* Validates that u, v, w views match MAC totals implied by the mesh and ghost width.
* Constructs a streamwise base profile (e.g. laminar or approximate turbulent shape) and applies random/structured perturbations controlled by `pert`.

Phase: **PreExchange**.

Future expansion: include the exact formulas used for TG and channel initial fields, and example configuration blocks.

---

## 6. Sub-grid-scale model (`SGS`)

`SGS` implements a **Smagorinsky** large-eddy model on the MAC grid:

* Parameters (from KV):

  * `dx`, `dy`, `dz` — cell sizes (typically match mesh spacing).
  * `Cs` — Smagorinsky constant (default ~0.17).
* Requires u, v, w to be registered and MAC-staggered.
* Optionally writes a cell-centred eddy viscosity field `nu_t` if present in the `FieldCatalog`.
* Uses `numerics::kernels::sgs_smagorinsky` to compute `nu_t` at cell centres from MAC velocities.

Phase: **Interior** — runs on interior cells after halos are current.

This part of the program and the kernel it references is where the heart of SGS model research is meant to be conducted.

---

## 7. Momentum predictor (`MomentumPredictor`)

`MomentumPredictor` advances the MAC velocity **without** imposing incompressibility (i.e. it builds a provisional velocity `u*`).

### 7.1 Time-integration modes

The predictor supports multiple temporal schemes, selected via `time_scheme`:

* `"fe"` — Forward Euler (explicit).
* `"be"` / `"backward_euler"` — first-order implicit.
* `"abm3"`, `"ab3"`, `"ab3am3"` — third-order Adams–Bashforth–Moulton predictor–corrector.

It keeps short histories of advection and diffusion terms to support multi-step methods.

### 7.2 Advection & diffusion

* Nonlinear advection is computed on the MAC grid using the KK3 scheme from the numerics library.
* Diffusion can be:

  * Explicit (Forward Euler diffusive step), or
  * Implicit using backward Euler with inner iterations.

Implicit controls:

* `pred_imp_max_iters` — maximum inner iterations.
* `pred_imp_rtol` — relative tolerance.
* `pred_imp_solver` — inner solver (`"jacobi"` or `"rbgs"`).
* `pred_imp_order` — order of the implicit stencil variant (1–3).

### 7.3 Forcing and options

Further controls include:

* `advect` — `on/off` toggle for the advection term.
* `fx`, `fy`, `fz` — constant body-force components.

The predictor:

* Validates that MAC field extents match what the mesh + ghost width imply.
* Requires at least two ghost cells (`ng >= 2`) when advection is enabled (KK3 needs a wider stencil).

Phase: **Interior**.

---

## 8. Pressure Poisson solve (`PressurePoisson`)

The pressure step enforces incompressibility by solving a **pressure Poisson equation** built from the divergence of `u*`.

Key design points:

* Uses PETSc (DMDA + KSP + PCMG) in a **matrix-free** configuration for the fine-level operator (`Amat`).
* Assembles a fine-grid AIJ matrix for `Pmat`, and lets PCMG build coarse levels via **Galerkin** coarsening.
* Uses cell-centred DMDA (`DMDA_Q0`) and interpolation for multigrid transfers.
* Supports both **constant-density** and **variable-density** cases via `β = 1/ρ` stored on cells.

Boundary conditions:

* Neumann (natural flux) and Dirichlet (fixed pressure) are supported.
* Periodic directions are deduced from BC specs and/or mesh periodicity.
* `mirror` on pressure is treated like homogeneous Neumann.

Operators:

* `ShellMult` implements the matrix-free stencil for `y = A x`, using face transmissibilities and cell volumes.
* A cached diagonal is built for use by Jacobi/Chebyshev smoothers.

Solver:

* Multigrid preconditioner with Jacobi–Chebyshev smoothing on each level.
* Outer KSP scheme:

  * PiPeCG for well-posed problems.
  * MINRES for singular configurations (e.g. pure Neumann).

This section can be expanded with:

* Exact discrete operator formula.
* Config keys (KSP type, tolerances, maximum iterations, monitoring options).
* Notes on singular cases and pressure-nullspace handling.

---

## 9. Velocity corrector (`VelocityCorrector`)

After solving for pressure, `VelocityCorrector` applies the **projection step**:

* Requires `u`, `v`, `w`, and `p` to be registered.
* Uses numerics kernels to:

  1. Compute face-centred pressure gradients from cell-centred pressure.
  2. Update MAC velocities in place using either constant or variable density.
* Validates that MAC field extents match mesh-implied totals.

Parameters:

* `rho` — density (for constant-density projection).
* `dx`, `dy`, `dz` — spacings, used in the gradient.

Phase: **PostBC** — runs after BCs, ensuring ghosts are valid before the final corrected velocity is used elsewhere.

---

## 10. Projection loop (`ProjectionLoop`)

`ProjectionLoop` is a composite action that orchestrates one full **projection timestep**:

High-level structure (conceptual):

1. (Optional) SGS update.
2. Apply BCs to u, v, w, p.
3. Perform halo exchanges.
4. Run the momentum predictor to obtain provisional velocity `u*`.
5. Construct the divergence field and pressure RHS.
6. Solve the pressure Poisson equation.
7. Correct velocities using the new pressure.
8. Optionally compute and log divergence norms (e.g. L∞ of `∇·u`).

It owns its own `Options` object describing:

* Spatial steps `dx`, `dy`, `dz`.
* Tolerances for divergence and pressure solves.

Phase: **PostExchange** — assumes faces and ghosts are up-to-date from the predictor and BC steps.

---

## 11. Build & optimisation (CMake)

From `CMakeLists.txt`:

* Target: `physics_fluids`.
* Languages: C++ and Fortran (Fortran via numerics and PETSc if needed).
* Links against:

  * `core`
  * `numerics`
  * `PETSC::petsc`
  * `OpenMP::OpenMP_CXX`
* Compile options:

  * `-O3 -march=native` for C++.
* Link options:

  * On Linux, link with `-Wl,-z,defs` to catch undefined symbols at link time.
* RPATHs are set so that the plugin can locate PETSc and sibling libraries at runtime.

This section can be extended with notes on recommended compilers, PETSc build options, and debugging flags.

---

## 12. Configuration & extension points

### 12.1 Common configuration keys

This plugin is driven almost entirely by key–value configuration (e.g. from a TOML/INI file). Common keys include:

* Geometry and fluid:

  * `dx`, `dy`, `dz`
  * `rho` (density)
  * `nu` (molecular viscosity)
  * `Cs` (Smagorinsky constant)
* Time stepping & stability:

  * `time_scheme` (`fe`, `be`, `abm3`, ...)
  * Predictor implicit controls (`pred_imp_*`).
* Advection & forcing:

  * `advect = on|off`
  * `fx`, `fy`, `fz` (body forces)
* Boundary conditions:

  * `bc.u.north = dirichlet:0.0`
  * `bc.v.* = mirror:0.0`
  * `bc.p.west = neumann:0.0`

### 12.2 Adding new physics

To extend the fluids plugin:

1. Implement a new `IAction` in `include/fluids/` + `src/` (for example, a new turbulence model or source term).
2. Wire it into `Register.cpp` with a descriptive name.
3. Optionally, update `Program.cpp` to construct your action as part of the standard pipeline.
4. Document any new configuration keys here.
