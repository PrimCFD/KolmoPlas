#!/usr/bin/env python3
"""
Lid-driven cavity regression vs Ghia et al. (1982) for KolmoPlas.

Features
--------
* Reads 2D KolmoPlas cavity result from CGNS or HDF5:
    - CGNS: u_cell, v_cell from last *_Cell FlowSolution_t.
    - HDF5: u_cell, v_cell from last Step_xxxxxx group.
* Extracts:
    - u-velocity along vertical centerline x = 0.5
    - v-velocity along horizontal centerline y = 0.5
* Loads Ghia et al. benchmark data (TABLE I and II) from simple ASCII files
  compatible with public ghiau.txt / ghiav.txt gists.
* Interpolates your profiles to Ghia sample points, prints L2/L∞ errors,
  and overlays plots.
* Saves figures and data next to the result file by default.

Usage
-----
  python ldc_ghia.py \
      [--file cavity.cgns|cavity.h5] \
      [--data-dir DIR] \
      [--Re 1000] \
      [--ghia-u ghiau.txt] \
      [--ghia-v ghiav.txt] \
      [--out-prefix PREFIX] \
      [--formats png,svg,pdf] \
      [--no-plot]
"""

import argparse
from pathlib import Path
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt

import CGNS.MAP
import CGNS.PAT.cgnskeywords as CK

try:
    import h5py
except ImportError:
    h5py = None


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS_DIR = REPO_ROOT / "build-regression" / "tests" /  "regression" / "cavity256" / "out" # adjust as needed, e.g. REPO_ROOT / "build"

DEFAULT_GHIA_DIR = REPO_ROOT / "scripts" / "ghia"

DEFAULT_GHIA_U = DEFAULT_GHIA_DIR / "ghiau.txt"
DEFAULT_GHIA_V = DEFAULT_GHIA_DIR / "ghiav.txt"

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

# ---------------------------------------------------------------------------
# CGNS helpers (aligned with tgv_spectrum.py style)
# ---------------------------------------------------------------------------

def normalize_label(x):
    """Convert CGNS label (numpy array / bytes / str) to a plain string."""
    import numpy as _np
    if isinstance(x, _np.ndarray):
        return x.tobytes().decode("ascii", errors="ignore").strip("\x00")
    if isinstance(x, (bytes, bytearray)):
        return x.decode("ascii", errors="ignore").strip("\x00")
    return str(x)


def find_all_children_of_type(node, cgns_type):
    target = normalize_label(cgns_type)
    return [c for c in node[2] if normalize_label(c[3]) == target]


def find_dataarray(node, name):
    for child in node[2]:
        if (normalize_label(child[3]) == normalize_label(CK.DataArray_ts)
                and child[0] == name):
            return child
    raise KeyError(f"DataArray_t '{name}' not found under node {node[0]!r}")


def pick_last_cell_solution(zone_node):
    sols = find_all_children_of_type(zone_node, CK.FlowSolution_ts)
    if not sols:
        raise RuntimeError("No FlowSolution_t nodes in zone.")
    cell_sols = [s for s in sols if s[0].endswith("_Cell")]
    if not cell_sols:
        return sorted(sols, key=lambda n: n[0])[-1]
    return sorted(cell_sols, key=lambda n: n[0])[-1]


def load_uv_from_cgns(filename: str) -> Tuple[np.ndarray, np.ndarray]:
    tree, links, paths = CGNS.MAP.load(filename)

    bases = find_all_children_of_type(tree, CK.CGNSBase_ts)
    if not bases:
        raise RuntimeError("No CGNSBase_t found.")
    base = bases[0]

    zones = find_all_children_of_type(base, CK.Zone_ts)
    if not zones:
        raise RuntimeError("No Zone_t found.")
    zone = zones[0]

    sol = pick_last_cell_solution(zone)
    print(f"[ldc_ghia] Using FlowSolution: {sol[0]}")

    u_node = find_dataarray(sol, "u_cell")
    v_node = find_dataarray(sol, "v_cell")

    # Data arrays live in node[1]
    u = np.array(u_node[1], dtype=np.float64)
    v = np.array(v_node[1], dtype=np.float64)

    if u.shape != v.shape:
        raise RuntimeError(f"Shape mismatch in CGNS u_cell/v_cell: {u.shape} vs {v.shape}")

    # Allow 2D or 3D (Nx, Ny, [Nz])
    return u, v


# ---------------------------------------------------------------------------
# HDF5 helpers (XdmfHdf5Writer layout)
# ---------------------------------------------------------------------------

def load_uv_from_hdf5(filename: str) -> Tuple[np.ndarray, np.ndarray]:
    if h5py is None:
        raise RuntimeError("h5py not available; cannot read HDF5 files.")

    with h5py.File(filename, "r") as f:
        step_names = sorted([k for k in f.keys() if k.startswith("Step_")])
        if not step_names:
            raise RuntimeError("No 'Step_xxxxxx' groups found in HDF5 file.")
        last_step = step_names[-1]
        g = f[last_step]
        if "u_cell" not in g or "v_cell" not in g:
            raise RuntimeError(f"Step group '{last_step}' missing u_cell/v_cell datasets.")

        u = np.array(g["u_cell"], dtype=np.float64)
        v = np.array(g["v_cell"], dtype=np.float64)

        if u.shape != v.shape:
            raise RuntimeError(f"Shape mismatch in HDF5 u_cell/v_cell: {u.shape} vs {v.shape}")

    print(f"[ldc_ghia] Using HDF5 step: {last_step}")
    return u, v


# ---------------------------------------------------------------------------
# Ghia data loaders
# ---------------------------------------------------------------------------

def parse_ghia_table(path: Path, Re: int) -> Tuple[np.ndarray, np.ndarray]:
    lines = path.read_text().splitlines()
    header_cols = None

    # --- Find the actual table header line ---
    for line in lines:
        line = line.strip()
        if not line.startswith("#"):
            continue
        # Strip '#' and split; we want lines where the first *token* is 'y' or 'x'
        tokens = line.lstrip("#").split()
        if not tokens:
            continue
        if tokens[0] in ("y", "x"):
            # e.g. "#  y      100      400      1000 ..."
            header_cols = tokens
            break

    if header_cols is None:
        raise RuntimeError(f"Could not find header line with Re columns in {path}")

    # First entry is "y" or "x"; the rest are Reynolds numbers.
    try:
        col_index = header_cols.index(str(Re))
    except ValueError:
        raise RuntimeError(
            f"Requested Re={Re} not found in header {header_cols} of {path}"
        )

    coord = []
    val = []

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) <= col_index:
            continue
        coord.append(float(parts[0]))
        val.append(float(parts[col_index]))

    return np.array(coord, dtype=float), np.array(val, dtype=float)



def load_ghia_u(path: Path, Re: int) -> Tuple[np.ndarray, np.ndarray]:
    return parse_ghia_table(path, Re)


def load_ghia_v(path: Path, Re: int) -> Tuple[np.ndarray, np.ndarray]:
    return parse_ghia_table(path, Re)


# ---------------------------------------------------------------------------
# Centerline extraction from solver field
# ---------------------------------------------------------------------------

def extract_centerlines(u: np.ndarray, v: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract:
      - u vs y along vertical centerline at x=0.5
      - v vs x along horizontal centerline at y=0.5

    Assumes:
      u, v shapes are (Nx, Ny, Nz) or (Nx, Ny) with unit box [0,1]x[0,1].
      If Nz>1, we take the mid-plane in z.
    """
    if u.ndim == 2:
        Nx, Ny = u.shape
        Nz = 1
        u2 = u[:, :, None]
        v2 = v[:, :, None]
    elif u.ndim == 3:
        Nx, Ny, Nz = u.shape
        u2 = u
        v2 = v
    else:
        raise RuntimeError(f"Unexpected velocity array dimension: {u.ndim}")

    i_mid = Nx // 2
    j_mid = Ny // 2
    k_mid = Nz // 2

    u_centerline = u2[i_mid, :, k_mid]   # u(x=0.5, y)
    v_centerline = v2[:, j_mid, k_mid]   # v(x, y=0.5)

    y_sim = (np.arange(Ny) + 0.5) / Ny
    x_sim = (np.arange(Nx) + 0.5) / Nx

    return x_sim, y_sim, v_centerline, u_centerline


# ---------------------------------------------------------------------------
# File discovery / dispatch (mirrors tgv_spectrum.py)
# ---------------------------------------------------------------------------

def guess_result_file(explicit_file: str, data_dir: Path) -> Path:
    if explicit_file:
        p = Path(explicit_file)
        if p.is_file():
            return p
        raise FileNotFoundError(f"Result file not found: {p}")

    candidates = []
    if data_dir.is_file():
        return data_dir

    if data_dir.is_dir():
        patterns = ("*.cgns", "*.cgns.h5", "*.h5", "*.hdf5")
        for pat in patterns:
            candidates.extend(sorted(data_dir.glob(pat)))

    if not candidates:
        raise FileNotFoundError(
            f"No result files (*.cgns, *.h5, *.hdf5) found under {data_dir!s}. "
            "Pass --file explicitly or adjust --data-dir."
        )

    latest = max(candidates, key=lambda p: p.stat().st_mtime)
    return latest


def load_u_v(filename: Path) -> Tuple[np.ndarray, np.ndarray]:
    ext = filename.suffix.lower()
    if ext == ".cgns":
        return load_uv_from_cgns(str(filename))
    elif ext in (".h5", ".hdf5"):
        return load_uv_from_hdf5(str(filename))
    else:
        raise RuntimeError(f"Unsupported file extension '{ext}' (expected .cgns/.h5/.hdf5)")


# ---------------------------------------------------------------------------
# Plotting, error metrics, and data export
# ---------------------------------------------------------------------------

def interp_profile(x_src: np.ndarray, y_src: np.ndarray, x_target: np.ndarray) -> np.ndarray:
    order = np.argsort(x_src)
    xs = x_src[order]
    ys = y_src[order]
    return np.interp(x_target, xs, ys)


def compute_errors(ref: np.ndarray, sim: np.ndarray):
    diff = sim - ref
    l2 = np.sqrt(np.mean(diff**2))
    linf = np.max(np.abs(diff))
    return l2, linf


def plot_centerlines(
    x_sim: np.ndarray,
    y_sim: np.ndarray,
    u_centerline: np.ndarray,
    v_centerline: np.ndarray,
    y_ghia: np.ndarray,
    u_ghia: np.ndarray,
    x_ghia: np.ndarray,
    v_ghia: np.ndarray,
    Re: int,
    label_prefix: str,
    out_prefix: str | None,
    formats: list[str],
    show: bool,
    u_l2: float,
    u_linf: float,
    v_l2: float,
    v_linf: float,
):

    # ---------- u(y) at x = 0.5 ----------
    fig_u, ax_u = plt.subplots(figsize=(5.0, 4.0))

    ax_u.plot(
        u_centerline, y_sim,
        "-o", label=label_prefix,
    )
    ax_u.plot(
        u_ghia, y_ghia,
        "s", mfc="none",
        label=f"Ghia et al. Re={Re}",
    )

    ax_u.set_xlabel(r"$u$")
    ax_u.set_ylabel(r"$y$")
    ax_u.set_title(rf"Lid-driven cavity: $u(y)$ at $x=0.5$, Re={Re}")
    ax_u.grid(True, ls="--", alpha=0.4)
    ax_u.legend(loc="best", frameon=False)

    # Add error annotation in the top-left corner
    text_u = (
        rf"$L_2$ = {u_l2:.3e}" + "\n" +
        rf"$L_\infty$ = {u_linf:.3e}"
    )
    ax_u.text(
        0.03, 0.97, text_u,
        transform=ax_u.transAxes,
        va="top", ha="left",
        bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.9),
    )

    fig_u.tight_layout()

    # ---------- v(x) at y = 0.5 ----------
    fig_v, ax_v = plt.subplots(figsize=(5.0, 4.0))

    ax_v.plot(
        x_sim, v_centerline,
        "-o", label=label_prefix,
    )
    ax_v.plot(
        x_ghia, v_ghia,
        "s", mfc="none",
        label=f"Ghia et al. Re={Re}",
    )

    ax_v.set_xlabel(r"$x$")
    ax_v.set_ylabel(r"$v$")
    ax_v.set_title(rf"Lid-driven cavity: $v(x)$ at $y=0.5$, Re={Re}")
    ax_v.grid(True, ls="--", alpha=0.4)
    ax_v.legend(loc="best", frameon=False)

    # Error annotation for v-profile
    text_v = (
        rf"$L_2$ = {v_l2:.3e}" + "\n" +
        rf"$L_\infty$ = {v_linf:.3e}"
    )
    ax_v.text(
        0.03, 0.97, text_v,
        transform=ax_v.transAxes,
        va="top", ha="left",
        bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.9),
    )

    fig_v.tight_layout()

    # ---------- Save ----------
    if out_prefix is not None:
        for ext in formats:
            ext = ext.strip().lower()
            if not ext:
                continue
            fig_u.savefig(f"{out_prefix}_u_profile.{ext}", bbox_inches="tight")
            fig_v.savefig(f"{out_prefix}_v_profile.{ext}", bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig_u)
        plt.close(fig_v)


def save_centerline_data(
    x_sim: np.ndarray,
    y_sim: np.ndarray,
    u_centerline: np.ndarray,
    v_centerline: np.ndarray,
    y_ghia: np.ndarray,
    u_ghia: np.ndarray,
    x_ghia: np.ndarray,
    v_ghia: np.ndarray,
    u_sim_on_ghia: np.ndarray,
    v_sim_on_ghia: np.ndarray,
    out_prefix: str,
):
    # u(y) data
    u_data = np.column_stack((y_sim, u_centerline))
    np.savetxt(
        f"{out_prefix}_u_profile_sim.dat",
        u_data,
        header="y_sim  u_sim",
        comments="",
    )

    u_table = np.column_stack((y_ghia, u_ghia, u_sim_on_ghia))
    np.savetxt(
        f"{out_prefix}_u_profile_ghia.dat",
        u_table,
        header="y_ghia  u_ghia  u_sim_on_ghia",
        comments="",
    )

    np.savetxt(
        f"{out_prefix}_u_profile.csv",
        np.column_stack((y_ghia, u_ghia, u_sim_on_ghia)),
        delimiter=",",
        header="y_ghia,u_ghia,u_sim_on_ghia",
        comments="",
    )

    # v(x) data
    v_data = np.column_stack((x_sim, v_centerline))
    np.savetxt(
        f"{out_prefix}_v_profile_sim.dat",
        v_data,
        header="x_sim  v_sim",
        comments="",
    )

    v_table = np.column_stack((x_ghia, v_ghia, v_sim_on_ghia))
    np.savetxt(
        f"{out_prefix}_v_profile_ghia.dat",
        v_table,
        header="x_ghia  v_ghia  v_sim_on_ghia",
        comments="",
    )

    np.savetxt(
        f"{out_prefix}_v_profile.csv",
        np.column_stack((x_ghia, v_ghia, v_sim_on_ghia)),
        delimiter=",",
        header="x_ghia,v_ghia,v_sim_on_ghia",
        comments="",
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Lid-driven cavity regression vs Ghia et al. (u(y), v(x))."
    )
    parser.add_argument(
        "--file",
        help="Explicit CGNS/HDF5 result file (overrides auto-discovery).",
        default=None,
    )
    parser.add_argument(
        "--data-dir",
        help=(
            "Directory to search for result files when --file is omitted. "
            f"Default: {DEFAULT_RESULTS_DIR}"
        ),
        default=None,
    )
    parser.add_argument(
        "--Re",
        type=int,
        default=1000,
        help="Reynolds number to compare against from Ghia tables (default: 1000).",
    )
    parser.add_argument(
        "--ghia-u",
        type=str,
        default=str(DEFAULT_GHIA_U),
        help=f"Ghia TABLE I file (ghiau.txt). Default: {DEFAULT_GHIA_U}",
    )
    parser.add_argument(
        "--ghia-v",
        type=str,
        default=str(DEFAULT_GHIA_V),
        help=f"Ghia TABLE II file (ghiav.txt). Default: {DEFAULT_GHIA_V}",
    )
    parser.add_argument(
        "--out-prefix",
        help=(
            "Output file prefix for saved plots and data. "
            "Default: result_file_stem + '_ghia_Re<Re>' in the result file directory."
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
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Compute errors and save plots/data only; do not show plots.",
    )

    args = parser.parse_args()

    data_dir = Path(args.data_dir) if args.data_dir else DEFAULT_RESULTS_DIR
    result_path = guess_result_file(args.file, data_dir)
    print(f"[ldc_ghia] Using result file: {result_path}")

    u, v = load_u_v(result_path)
    x_sim, y_sim, v_centerline, u_centerline = extract_centerlines(u, v)

    # Load Ghia reference data
    ghia_u_path = Path(args.ghia_u)
    ghia_v_path = Path(args.ghia_v)
    if not ghia_u_path.is_file():
        raise FileNotFoundError(f"Ghia u-table not found at {ghia_u_path}")
    if not ghia_v_path.is_file():
        raise FileNotFoundError(f"Ghia v-table not found at {ghia_v_path}")

    y_ghia, u_ghia = load_ghia_u(ghia_u_path, args.Re)
    x_ghia, v_ghia = load_ghia_v(ghia_v_path, args.Re)

    # Interpolate sim → Ghia points for error metrics
    u_sim_on_ghia = interp_profile(y_sim, u_centerline, y_ghia)
    v_sim_on_ghia = interp_profile(x_sim, v_centerline, x_ghia)

    u_l2, u_linf = compute_errors(u_ghia, u_sim_on_ghia)
    v_l2, v_linf = compute_errors(v_ghia, v_sim_on_ghia)

    print(f"[ldc_ghia] Re = {args.Re}")
    print("  u(y) centerline errors:  L2 = {:.3e}, Linf = {:.3e}".format(u_l2, u_linf))
    print("  v(x) centerline errors:  L2 = {:.3e}, Linf = {:.3e}".format(v_l2, v_linf))

    # Determine output prefix (for plots + data)
    if args.out_prefix is not None:
        out_prefix = args.out_prefix
    else:
        stem = result_path.with_suffix("").name
        out_prefix = str(result_path.parent / f"{stem}_ghia_Re{args.Re}")

    # Save data tables
    save_centerline_data(
        x_sim,
        y_sim,
        u_centerline,
        v_centerline,
        y_ghia,
        u_ghia,
        x_ghia,
        v_ghia,
        u_sim_on_ghia,
        v_sim_on_ghia,
        out_prefix,
    )

    # Parse formats
    formats = [s.strip() for s in args.formats.split(",") if s.strip()]

    # Make / save plots
    show_plots = not args.no_plot
    label_prefix = f"KolmoPlas ({result_path.name})"
    plot_centerlines(
        x_sim,
        y_sim,
        u_centerline,
        v_centerline,
        y_ghia,
        u_ghia,
        x_ghia,
        v_ghia,
        args.Re,
        label_prefix,
        out_prefix,
        formats,
        show_plots,
        u_l2,
        u_linf,
        v_l2,
        v_linf,
    )


if __name__ == "__main__":
    main()
