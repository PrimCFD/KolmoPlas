#!/usr/bin/env python3
"""
Compute an isotropic energy spectrum E(k) from a CGNS or HDF5 result file
written by this solver's CGNSWriter/XdmfHdf5Writer (cell-centered u_cell, v_cell, w_cell).

Usage:
    python tgv_spectrum.py result.cgns [--solution-name NAME] [--no-plot]

Usage:
    python tgv_spectrum.py \
        [--file result.cgns|result.h5] \
        [--data-dir DIR] \
        [--solution-name NAME] \
        [--step-index N] \
        [--no-plot]

By default, the script:
  * for CGNS: picks the last FlowSolutionAtStepXXXX_Cell
  * for HDF5: picks the last Step_xxxxxx group
  * expects u_cell, v_cell, w_cell at CellCenter
  * assumes a periodic uniform box with Lx=NX, Ly=NY, Lz=NZ.
"""

import argparse
from pathlib import Path
import os

import numpy as np
import matplotlib.pyplot as plt

import CGNS.MAP
import CGNS.PAT.cgnsutils as CGU
import CGNS.PAT.cgnskeywords as CK

try:
    import h5py
except ImportError:
    h5py = None

PAPER_STYLE = {
    "figure.figsize": (5.0, 4.0),
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.fontsize": 10,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "lines.linewidth": 1.0,
    "lines.markersize": 4,
}

plt.rcParams.update(PAPER_STYLE)

def _normalize_label(x):
    """
    Convert a CGNS label (which may be a numpy array, bytes, or str)
    into a plain Python string for reliable comparisons.
    """
    import numpy as _np

    # numpy array -> bytes -> string
    if isinstance(x, _np.ndarray):
        # pyCGNS usually stores labels as fixed-width byte arrays
        return x.tobytes().decode("ascii", errors="ignore").strip("\x00")

    # raw bytes -> string
    if isinstance(x, (bytes, bytearray)):
        return x.decode("ascii", errors="ignore").strip("\x00")

    # already a string or something printable
    return str(x)


def find_first_child_of_type(node, cgns_type):
    """
    node: CGNS/Python node [name, value, children, type]
    cgns_type: e.g. CK.CGNSBase_ts, CK.Zone_ts, CK.FlowSolution_ts

    Returns the first child node of that type.
    """
    target = _normalize_label(cgns_type)
    for child in node[2]:  # children list
        if _normalize_label(child[3]) == target:  # type is at index 3
            return child
    raise RuntimeError(f"No child of type {target!r} under node {node[0]!r}")


def find_all_children_of_type(node, cgns_type):
    """
    Return all children of given CGNS/SIDS type under `node`.
    """
    target = _normalize_label(cgns_type)
    return [c for c in node[2] if _normalize_label(c[3]) == target]



def find_dataarray(node, name):
    """
    Return the DataArray_t node with given name under `node`.

    CGNS/Python layout: [name, value, children, type]
    """
    for child in node[2]:
        if _normalize_label(child[3]) == _normalize_label(CK.DataArray_ts) and child[0] == name:
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

    # Data is stored in node[1] as a numpy array (Fortran order).
    u = np.array(u_node[1], dtype=np.float64)
    v = np.array(v_node[1], dtype=np.float64)
    w = np.array(w_node[1], dtype=np.float64)

    # Ensure 3D shape; pyCGNS usually gives correct shape already.
    if u.ndim != 3:
        raise RuntimeError(f"u_cell has unexpected shape {u.shape} (expected 3D).")

    return u, v, w


# ---------------------------------------------------------------------------
# HDF5 / multi-backend helpers
# ---------------------------------------------------------------------------

# Default directory used when no filename is given (can be overridden with --data-dir).
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS_DIR = REPO_ROOT / "build-regression" / "tests" /  "regression" / "tgv128" / "out" # adjust as needed, e.g. REPO_ROOT / "build"


def load_velocity_from_hdf5(filename, step_index=None):
    """
    Load u_cell, v_cell, w_cell from an HDF5 file written by XdmfHdf5Writer.

    Layout (see XdmfHdf5Writer.cpp):
      /Step_xxxxxx/u_cell, v_cell, w_cell  (or fallback to u, v, w)

    step_index:
      * None (default): use the last available Step_xxxxxx group.
      * >= 0: 0-based index into sorted Step_* groups.
      * < 0 : negative index from the end (-1 = last, etc.).
    """
    if h5py is None:
        raise RuntimeError("h5py is not available; cannot read HDF5 files.")

    with h5py.File(filename, "r") as f:
        step_names = sorted([name for name in f.keys() if name.startswith("Step_")])
        if not step_names:
            raise RuntimeError("No 'Step_xxxxxx' groups found in HDF5 file.")

        if step_index is None:
            idx = len(step_names) - 1
        elif step_index < 0:
            idx = len(step_names) + step_index
        else:
            idx = step_index

        if idx < 0 or idx >= len(step_names):
            raise IndexError(
                f"step_index {step_index} is out of range for {len(step_names)} steps."
            )

        step_name = step_names[idx]
        print(f"Using HDF5 group: {step_name}")
        g = f[step_name]

        def _get_component(alias, base):
            if alias in g:
                return np.array(g[alias], dtype=np.float64)
            if base in g:
                return np.array(g[base], dtype=np.float64)
            raise RuntimeError(
                f"Neither '{alias}' nor '{base}' found in HDF5 group '{step_name}'."
            )

        u = _get_component("u_cell", "u")
        v = _get_component("v_cell", "v")
        w = _get_component("w_cell", "w")

        if u.shape != v.shape or u.shape != w.shape:
            raise RuntimeError(
                f"Shape mismatch in HDF5 velocities: u{u.shape}, v{v.shape}, w{w.shape}"
            )

        if u.ndim != 3:
            raise RuntimeError(f"HDF5 velocity has unexpected shape {u.shape} (expected 3D).")

        return u, v, w


def guess_result_file(explicit_filename, data_dir: Path) -> Path:
    """
    Resolve the result file to use, with a small priority rule:

      1) If explicit_filename is provided and exists, use it.
      2) Otherwise, search data_dir for *.cgns, *.cgns.h5, *.h5, *.hdf5
         and pick the most recently modified file.

    Raises FileNotFoundError if nothing is found.
    """
    if explicit_filename:
        p = Path(explicit_filename)
        if p.is_file():
            return p
        raise FileNotFoundError(f"Result file not found: {p}")

    data_dir = Path(data_dir)
    if data_dir.is_file():
        return data_dir

    candidates = []
    if data_dir.is_dir():
        patterns = ("*.cgns", "*.cgns.h5", "*.h5", "*.hdf5")
        for pat in patterns:
            candidates.extend(sorted(data_dir.glob(pat)))

    if not candidates:
        raise FileNotFoundError(
            f"No result files (*.cgns, *.h5, *.hdf5) found under {data_dir!s}."
        )

    # Pick the newest candidate.
    latest = max(candidates, key=lambda p: p.stat().st_mtime)
    return latest


def load_velocity(filename: Path, solution_name=None, step_index=None):
    """
    Dispatch to the appropriate backend based on filename extension.
    """
    suffix = filename.suffix.lower()
    if suffix == ".cgns":
        return load_velocity_from_cgns(str(filename), solution_name=solution_name)
    elif suffix in (".h5", ".hdf5"):
        return load_velocity_from_hdf5(str(filename), step_index=step_index)
    else:
        raise RuntimeError(
            f"Unsupported file extension '{suffix}' (expected .cgns, .h5, or .hdf5)."
        )


def compute_isotropic_spectrum(u, v, w, nbins=64, log_bins=True, min_modes=0):
    """
    Compute isotropic 1D spectrum E(k) from 3D velocity fields u,v,w (cell-centered).

    Steps:
      1) subtract mean
      2) 3D FFT for each component
      3) kinetic energy density Ef(kx,ky,kz) = 0.5*(|û|^2 + |v̂|^2 + |ŵ|^2)
      4) bin Ef into spherical shells in k-space

    Parameters
    ----------
    nbins : int
        Number of k-shells.
    log_bins : bool
        If True (default), use logarithmic spacing in k.
        If False, use linear spacing.
    min_modes : int
        Minimum number of Fourier modes contributing to a bin.
        Bins with fewer contributing modes are dropped.
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

    # Flatten arrays
    k_flat = k_mag.ravel()
    E_flat = E_kxyz.ravel()

    # Only positive k for spectrum
    mask_pos = k_flat > 0.0
    k_pos = k_flat[mask_pos]
    E_pos = E_flat[mask_pos]

    k_max = k_pos.max()
    k_min = k_pos.min()

    if nbins is None or nbins <= 0:
        nbins = 64

    if log_bins:
        # Logarithmic bins: edges from k_min to k_max
        bin_edges = np.logspace(np.log10(k_min), np.log10(k_max), nbins + 1)
        # Geometric mean for bin center
        bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])
    else:
        # Linear bins
        bin_edges = np.linspace(k_min, k_max, nbins + 1)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Bin indices: 0 .. nbins-1
    bin_idx = np.digitize(k_pos, bin_edges) - 1
    # Clamp any out-of-range (numerical noise) to last bin
    bin_idx = np.clip(bin_idx, 0, nbins - 1)

    shell_energy = np.bincount(bin_idx, weights=E_pos, minlength=nbins)
    shell_count = np.bincount(bin_idx, minlength=nbins)

    # First remove bins that are completely empty
    nonempty = shell_count > 0
    k_vals = bin_centers[nonempty]
    E1D = shell_energy[nonempty] / shell_count[nonempty]
    modes = shell_count[nonempty]

    # Drop obviously empty very-low-energy bins AND bins with too few modes
    E_max = E1D.max()
    tiny_thresh = 1e-12 * E_max  # only truly negligible bins
    min_modes = max(int(min_modes), 1)

    keep = (E1D > tiny_thresh) & (modes >= min_modes)

    k_vals = k_vals[keep]
    E1D = E1D[keep]

    return k_vals, E1D

def plot_spectrum(k_vals, E1D, out_prefix=None, formats=None, show=True):
    """
    Make 'paper-ready' basic and compensated spectra plots.

    If out_prefix is not None, save each figure as:
      <out_prefix>_spectrum.<ext>
      <out_prefix>_spectrum_compensated.<ext>
    for each ext in `formats` (list of strings).

    `show` controls whether plt.show() is called.
    """
    if formats is None:
        formats = []

    # ------------------------------------------------------------------
    # 1) Basic spectrum: E(k) with fitted -5/3 reference
    # ------------------------------------------------------------------
    fig1, ax1 = plt.subplots()

    # Simulation data as points only
    ax1.loglog(
        k_vals,
        E1D,
        marker="o",
        linestyle="None",
        label="simulation",
    )

    # Fit an inertial-range k^{-5/3} line in the central part of the spectrum
    if len(k_vals) >= 10:
        s = -5.0 / 3.0

        # Fit over the central 50% of available k (avoid extremes)
        i0 = int(len(k_vals) * 0.25)
        i1 = int(len(k_vals) * 0.75)
        i0 = max(0, min(i0, len(k_vals) - 2))
        i1 = max(i0 + 1, min(i1, len(k_vals) - 1))

        k_fit = k_vals[i0:i1]
        E_fit = E1D[i0:i1]

        logk = np.log(k_fit)
        logE = np.log(E_fit)

        # log E ≈ log C + s log k  ⇒  log C = mean(logE − s logk)
        logC = np.mean(logE - s * logk)
        C = np.exp(logC)

        k_line = np.logspace(np.log10(k_fit[0]), np.log10(k_fit[-1]), 100)
        E_line = C * k_line**s

        ax1.loglog(
            k_line,
            E_line,
            linestyle="--",
            color="0.3",
            label=r"$k^{-5/3}$ fit",
        )

    # Axis labels / title
    ax1.set_xlabel(r"$k$")
    ax1.set_ylabel(r"$E(k)$")
    ax1.set_title("Isotropic energy spectrum")

    ax1.grid(True, which="both", linestyle="--", alpha=0.4)
    ax1.legend(loc="upper right", frameon=False)

    fig1.tight_layout()

    # ------------------------------------------------------------------
    # 2) Compensated spectrum: E(k) k^{5/3}
    # ------------------------------------------------------------------
    fig2, ax2 = plt.subplots()
    Ek_comp = E1D * (k_vals ** (5.0 / 3.0))

    ax2.semilogx(
        k_vals,
        Ek_comp,
        marker="o",
        linestyle="None",
    )

    ax2.set_xlabel(r"$k$")
    ax2.set_ylabel(r"$E(k)\,k^{5/3}$")
    ax2.set_title("Compensated spectrum")

    ax2.grid(True, which="both", linestyle="--", alpha=0.4)
    fig2.tight_layout()

    # ------------------------------------------------------------------
    # 3) Save figures if requested
    # ------------------------------------------------------------------
    if out_prefix is not None:
        for ext in formats:
            ext = ext.strip().lower()
            if not ext:
                continue
            fig1.savefig(f"{out_prefix}_spectrum.{ext}", bbox_inches="tight")
            fig2.savefig(f"{out_prefix}_spectrum_compensated.{ext}", bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig1)
        plt.close(fig2)



def save_spectrum_data(k_vals, E1D, out_prefix):
    """
    Save spectrum data for post-processing / papers.

    Creates:
      <out_prefix>_spectrum.csv  (k,E(k),Ek*k^(5/3) with CSV commas)
      <out_prefix>_spectrum.dat  (ASCII space-separated, columns: k E(k) Ek*k^(5/3))
    """
    Ek_comp = E1D * (k_vals ** (5.0 / 3.0))
    data = np.column_stack((k_vals, E1D, Ek_comp))

    # CSV (good for spreadsheets, etc.)
    csv_name = f"{out_prefix}_spectrum.csv"
    np.savetxt(csv_name, data, delimiter=",", header="k,E(k),E(k)*k^(5/3)", comments="")

    # Plain ASCII (space-separated), nice for gnuplot, etc.
    dat_name = f"{out_prefix}_spectrum.dat"
    np.savetxt(dat_name, data, header="k  E(k)  E(k)*k^(5/3)", comments="")



def main():
    parser = argparse.ArgumentParser(
        description="Kolmogorov-like energy spectrum from Taylor–Green result (CGNS or HDF5)."
    )
    parser.add_argument(
        "filename",
        nargs="?",
        help="Result file (CGNS or HDF5). Deprecated; prefer --file.",
        default=None,
    )
    parser.add_argument(
        "--file",
        dest="file",
        help=(
            "Explicit CGNS/HDF5 result file. "
            "If omitted, auto-discover the newest file in --data-dir."
        ),
        default=None,
    )
    parser.add_argument(
        "--data-dir",
        help=(
            "Directory to search for result files when --file/filename is omitted. "
            f"Default: {DEFAULT_RESULTS_DIR}"
        ),
        default=str(DEFAULT_RESULTS_DIR),
    )
    parser.add_argument(
        "--solution-name",
        help=(
            "For CGNS: exact FlowSolution_t name to use "
            "(default: last *'_Cell' solution, e.g. FlowSolutionAtStep000100_Cell)."
        ),
        default=None,
    )
    parser.add_argument(
        "--step-index",
        type=int,
        default=None,
        help=(
            "For HDF5: 0-based index into sorted Step_xxxxxx groups. "
            "Negative values count from the end (-1 = last). "
            "Default: last available step."
        ),
    )

    parser.add_argument(
        "--num-bins",
        type=int,
        default=64,
        help=(
            "Number of spherical k-shell bins for the isotropic spectrum "
            "(default: 64)."
        ),
    )

    parser.add_argument(
        "--linear-bins",
        help="Use linearly spaced k-shells instead of log-spaced.",
        action="store_true",
    )

    parser.add_argument(
        "--min-modes-per-bin",
        type=int,
        default=20,
        help=(
            "Minimum number of Fourier modes contributing to a k-bin; "
            "bins with fewer modes are discarded (default: 20)."
        ),
    )

    parser.add_argument(
        "--no-plot",
        help="Compute and print spectrum but do not pop up plots.",
        action="store_true",
    )

    parser.add_argument(
        "--out-prefix",
        help=(
            "Output file prefix for saved plots and data. "
            "Default: derived from result filename (e.g. tgv_128_spectrum)."
        ),
        default=None,
    )
    parser.add_argument(
        "--formats",
        help=(
            "Comma-separated list of image formats to save "
            "(e.g. 'png,svg,pdf'). Default: 'png,svg,pdf'."
        ),
        default="png,svg,pdf",
    )

    args = parser.parse_args()

    # Resolve result file
    explicit = args.file or args.filename
    data_dir = Path(args.data_dir)
    result_path = guess_result_file(explicit, data_dir)
    print(f"Using result file: {result_path}")

    # Load velocity field and compute spectrum
    u, v, w = load_velocity(result_path, solution_name=args.solution_name, step_index=args.step_index)
    k_vals, E1D = compute_isotropic_spectrum(u, v, w, nbins=args.num_bins, log_bins=not args.linear_bins, min_modes=args.min_modes_per_bin)

    # Print a small table of k and E(k)
    print("# k   E(k)")
    for k, e in zip(k_vals[:20], E1D[:20]):  # first 20 modes
        print(f"{k:10.5f}  {e:15.8e}")

    # Determine output prefix (for plots + data)
    if args.out_prefix is not None:
        # User-specified prefix: respect as-is (relative to CWD or absolute)
        out_prefix = args.out_prefix
    else:
        # Default: same directory as the result file, plus "_spectrum" suffix
        stem = result_path.with_suffix("").name
        out_prefix = str(result_path.parent / f"{stem}_spectrum")


    # Parse formats
    formats = [s.strip() for s in args.formats.split(",") if s.strip()]

    # Save spectrum data
    save_spectrum_data(k_vals, E1D, out_prefix)

    # Make / save plots
    show_plots = not args.no_plot
    plot_spectrum(k_vals, E1D, out_prefix=out_prefix, formats=formats, show=show_plots)



if __name__ == "__main__":
    main()
