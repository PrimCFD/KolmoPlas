#!/usr/bin/env python3
"""
Compute an isotropic energy spectrum E(k) from a CGNS result file
written by this solver's CGNSWriter (cell-centered u_cell, v_cell, w_cell).

Usage:
    python tgv_spectrum.py result.cgns [--solution-name NAME] [--no-plot]

By default, the script:
  * picks the last FlowSolutionAtStepXXXX_Cell
  * expects u_cell, v_cell, w_cell at CellCenter
  * assumes a periodic uniform box with Lx=NX, Ly=NY, Lz=NZ.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

import CGNS.MAP
import CGNS.PAT.cgnsutils as CGU
import CGNS.PAT.cgnskeywords as CK


def find_first_child_of_type(node, cgns_type):
    """
    node: CGNS/Python node [name, label, children, data]
    cgns_type: e.g. CK.CGNSBase_ts, CK.Zone_ts, CK.FlowSolution_ts

    Returns the first child node of that type.
    """
    for child in node[2]:
        if child[1] == cgns_type:
            return child
    raise RuntimeError(f"No child of type {cgns_type} under node {node[0]!r}")


def find_all_children_of_type(node, cgns_type):
    return [c for c in node[2] if c[1] == cgns_type]


def find_dataarray(node, name):
    """
    Return the DataArray_t node with given name under `node`.
    """
    for child in node[2]:
        if child[1] == CK.DataArray_ts and child[0] == name:
            return child
    raise KeyError(f"DataArray_t '{name}' not found under node {node[0]!r}")


def pick_solution_node(zone_node, explicit_name=None):
    """
    From a Zone_t node, pick a FlowSolution_t:
      - if explicit_name is given, use that (exact match)
      - otherwise, pick the last FlowSolution whose name endswith '_Cell'
        (FlowSolutionAtStep%06d_Cell from your writer).
    """
    sols = find_all_children_of_type(zone_node, CK.FlowSolution_ts)
    if not sols:
        raise RuntimeError("No FlowSolution_t nodes found in zone.")

    if explicit_name is not None:
        for s in sols:
            if s[0] == explicit_name:
                return s
        raise RuntimeError(f"FlowSolution_t named '{explicit_name}' not found.")

    # Auto mode: last *_Cell solution by name
    cell_sols = [s for s in sols if s[0].endswith("_Cell")]
    if not cell_sols:
        # fallback: just the last solution
        return sorted(sols, key=lambda n: n[0])[-1]

    return sorted(cell_sols, key=lambda n: n[0])[-1]


def load_velocity_from_cgns(filename, solution_name=None):
    """
    Load u_cell, v_cell, w_cell from the given CGNS file,
    returning u, v, w as 3D numpy arrays (double).
    """
    tree, links, paths = CGNS.MAP.load(filename)

    # Base + Zone (assume one each, which matches your writer)
    bases = find_all_children_of_type(tree, CK.CGNSBase_ts)
    if not bases:
        raise RuntimeError("No CGNSBase_t found in tree.")
    base = bases[0]

    zones = find_all_children_of_type(base, CK.Zone_ts)
    if not zones:
        raise RuntimeError("No Zone_t found in base.")
    zone = zones[0]

    sol = pick_solution_node(zone, explicit_name=solution_name)
    print(f"Using FlowSolution node: {sol[0]}")

    # DataArray_t nodes for cell-centered velocity aliases
    u_node = find_dataarray(sol, "u_cell")
    v_node = find_dataarray(sol, "v_cell")
    w_node = find_dataarray(sol, "w_cell")

    # Data is stored in node[3] as a numpy array (Fortran order).
    u = np.array(u_node[3], dtype=np.float64, copy=False)
    v = np.array(v_node[3], dtype=np.float64, copy=False)
    w = np.array(w_node[3], dtype=np.float64, copy=False)

    # Ensure 3D shape; pyCGNS usually gives correct shape already.
    if u.ndim != 3:
        raise RuntimeError(f"u_cell has unexpected shape {u.shape} (expected 3D).")

    return u, v, w


def compute_isotropic_spectrum(u, v, w):
    """
    Compute isotropic 1D spectrum E(k) from 3D velocity fields u,v,w (cell-centered).
    Algorithm:
      1) subtract mean
      2) 3D FFT for each component
      3) kinetic energy density Ef(kx,ky,kz) = 0.5*(|û|^2 + |v̂|^2 + |ŵ|^2)
      4) bin Ef into spherical shells in k-space (using integer |k| index)
    """
    # Remove means (should be ~0 for TGV, but safe)
    u_prime = u - u.mean()
    v_prime = v - v.mean()
    w_prime = w - w.mean()

    Nx, Ny, Nz = u_prime.shape
    print(f"Grid size: {Nx} x {Ny} x {Nz}")

    # 3D FFTs
    u_hat = np.fft.fftn(u_prime)
    v_hat = np.fft.fftn(v_prime)
    w_hat = np.fft.fftn(w_prime)

    # Spectral kinetic energy density (unnormalized; overall constant not important for slope)
    E_kxyz = 0.5 * (
        np.abs(u_hat) ** 2 +
        np.abs(v_hat) ** 2 +
        np.abs(w_hat) ** 2
    )

    # Wavenumber components (assume unit spacing, periodic box L=N)
    kx = 2.0 * np.pi * np.fft.fftfreq(Nx, d=1.0)
    ky = 2.0 * np.pi * np.fft.fftfreq(Ny, d=1.0)
    kz = 2.0 * np.pi * np.fft.fftfreq(Nz, d=1.0)

    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    k_mag = np.sqrt(KX**2 + KY**2 + KZ**2)

    # Integer shell index based on |k|
    k_idx = np.rint(k_mag).astype(int)

    E_flat = E_kxyz.ravel()
    k_flat = k_idx.ravel()

    kmax = k_flat.max()
    shell_energy = np.bincount(k_flat, weights=E_flat, minlength=kmax + 1)
    shell_count = np.bincount(k_flat, minlength=kmax + 1)

    # Avoid division by zero
    valid = shell_count > 0
    k_vals = np.arange(kmax + 1, dtype=float)[valid]
    E1D = shell_energy[valid] / shell_count[valid]

    return k_vals, E1D


def plot_spectrum(k_vals, E1D):
    # Basic spectrum
    plt.figure()
    plt.loglog(k_vals, E1D, marker="o", linestyle="-")
    plt.xlabel(r"$k$")
    plt.ylabel(r"$E(k)$")
    plt.title("Isotropic energy spectrum")
    plt.grid(True, which="both", ls="--", alpha=0.5)

    # Compensated spectrum E(k) k^{5/3}
    plt.figure()
    Ek_comp = E1D * (k_vals ** (5.0 / 3.0))
    plt.semilogx(k_vals, Ek_comp, marker="o", linestyle="-")
    plt.xlabel(r"$k$")
    plt.ylabel(r"$E(k) k^{5/3}$")
    plt.title("Compensated spectrum")
    plt.grid(True, which="both", ls="--", alpha=0.5)

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Kolmogorov-like energy spectrum from CGNS Taylor–Green result."
    )
    parser.add_argument("filename", help="CGNS result file (from this solver's CGNSWriter)")
    parser.add_argument(
        "--solution-name",
        help=(
            "Exact FlowSolution_t name to use "
            "(default: last *'_Cell' solution, e.g. FlowSolutionAtStep000100_Cell)"
        ),
        default=None,
    )
    parser.add_argument(
        "--no-plot",
        help="Compute and print spectrum but do not pop up plots.",
        action="store_true",
    )

    args = parser.parse_args()

    u, v, w = load_velocity_from_cgns(args.filename, solution_name=args.solution_name)
    k_vals, E1D = compute_isotropic_spectrum(u, v, w)

    # Print a small table of k and E(k)
    print("# k   E(k)")
    for k, e in zip(k_vals[:20], E1D[:20]):  # first 20 modes
        print(f"{k:10.5f}  {e:15.8e}")

    if not args.no_plot:
        plot_spectrum(k_vals, E1D)


if __name__ == "__main__":
    main()
