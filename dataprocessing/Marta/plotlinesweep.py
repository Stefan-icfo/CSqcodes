# =============================== STITCHED LINE-SWEEP PLOT ===============================
# Loads runs 809, 810, 911 from a QCoDeS DB, standardizes units (x in mV, y in V),
# stitches them on a common grid, and plots one image:
#   x-axis: 740–820 mV
#   y-axis: -0.3 V to max(y) in run 911
# Overlaps are merged by max intensity.
# ----------------------------------------------------------------------------------------

import os
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ============================== USER SETTINGS ===========================================
DB_PATH = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v9.db"
OUT_DIR = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Figures"
RUN_IDS = [809, 810, 911]        # runs to stitch

# Plot limits
X_MIN_MV, X_MAX_MV = 740.0, 820.0
Y_MIN_V = -0.3
Y_MAX_FROM_RUN = 911             # choose which run defines y-maximum

# Grid resolution (increase for smoother image)
NX, NY = 801, 801                # ~0.1 mV in x

SAVE_NAME = f"stitched_linesweep_{RUN_IDS[0]}_{RUN_IDS[1]}_{RUN_IDS[2]}.png"
# ========================================================================================

os.makedirs(OUT_DIR, exist_ok=True)
qc.config["core"]["db_location"] = DB_PATH


# ==================================== HELPERS ===========================================
def _unit_to_mV(u: str) -> float:
    if not u: return 1.0
    s = u.strip().lower()
    if s in {"mv", "millivolt", "millivolts"}: return 1.0
    if s in {"v", "volt", "volts"}:           return 1e3
    if s in {"uv", "µv", "microvolt", "microvolts"}: return 1e-3
    return 1.0

def _unit_to_V(u: str) -> float:
    if not u: return 1.0
    s = u.strip().lower()
    if s in {"v", "volt", "volts"}:           return 1.0
    if s in {"mv", "millivolt", "millivolts"}: return 1e-3
    if s in {"uv", "µv", "microvolt", "microvolts"}: return 1e-6
    return 1.0

def load_2d_xy_z(run_id):
    """
    Robustly load a 2D sweep: returns x, y, z as 1D arrays plus meta dict.
    Works even if pdata[dep] lacks an 'axes' field.
    """
    ds = qc.load_by_id(run_id)  # or from qcodes.dataset import load_by_id
    pdata = ds.get_parameter_data()
    if not pdata:
        raise RuntimeError(f"Run {run_id}: empty parameter data.")

    dep = next(iter(pdata.keys()))
    block = pdata[dep]

    # --- infer axis names ---
    xname = yname = None
    ps = getattr(ds, "paramspecs", None)
    ps_dep = ps.get(dep) if ps else None

    # Preferred: ParamSpec.setpoints (tuple/list of names)
    if ps_dep is not None and getattr(ps_dep, "setpoints", None):
        sps = list(ps_dep.setpoints)
        if len(sps) >= 2:
            xname, yname = sps[0], sps[1]

    # Fallback: ParamSpec.depends_on (comma-separated)
    if (xname is None or yname is None) and ps_dep is not None and getattr(ps_dep, "depends_on", None):
        parts = [p.strip() for p in ps_dep.depends_on.split(",") if p.strip()]
        if len(parts) >= 2:
            xname, yname = parts[0], parts[1]

    # Last resort: first two non-dependent keys in the block
    if xname is None or yname is None:
        candidates = [k for k in block.keys() if k != dep and not k.endswith("_errors")]
        if len(candidates) < 2:
            raise RuntimeError(
                f"Run {run_id}: couldn't infer 2D axes. Block keys: {list(block.keys())}"
            )
        xname, yname = candidates[:2]

    # --- data arrays ---
    try:
        x = np.asarray(block[xname]).ravel()
        y = np.asarray(block[yname]).ravel()
        z = np.asarray(block[dep]).ravel()
    except KeyError as e:
        raise RuntimeError(
            f"Run {run_id}: expected keys not found. "
            f"Have {list(block.keys())}, wanted {xname}, {yname}, {dep}"
        ) from e

    if x.size == 0 or y.size == 0 or z.size == 0:
        raise RuntimeError(f"Run {run_id}: one of x/y/z is empty (sizes: {x.size}, {y.size}, {z.size}).")

    # --- labels/units ---
    def _get(pspec, name, attr, default):
        return getattr(pspec.get(name, None), attr, default) if pspec else default

    x_unit = _get(ps, xname, "unit", "")
    y_unit = _get(ps, yname, "unit", "")
    z_unit = _get(ps, dep,   "unit", "")

    x_label = _get(ps, xname, "label", xname)
    y_label = _get(ps, yname, "label", yname)
    z_label = _get(ps, dep,   "label", dep)

    meta = dict(
        x_name=xname, y_name=yname, z_name=dep,
        x_unit=x_unit, y_unit=y_unit, z_unit=z_unit,
        x_label=x_label, y_label=y_label, z_label=z_label,
        run_id=run_id, guid=getattr(ds, "guid", None),
    )
    return x, y, z, meta

def bin_to_grid(x_mV, y_V, z, x_edges, y_edges, reducer="max"):
    """Bin scattered (x,y,z) to a regular grid. reducer='max' or 'mean'."""
    xi = np.digitize(x_mV, x_edges) - 1
    yi = np.digitize(y_V,  y_edges) - 1
    nxg, nyg = len(x_edges) - 1, len(y_edges) - 1

    mask = (xi >= 0) & (xi < nxg) & (yi >= 0) & (yi < nyg) & np.isfinite(z)
    if not np.any(mask):
        return np.full((nyg, nxg), np.nan)

    xi, yi, z = xi[mask], yi[mask], z[mask]
    grid = np.full((nyg, nxg), np.nan)

    if reducer == "max":
        lin = yi * nxg + xi
        order = np.argsort(lin, kind="mergesort")
        lin, z = lin[order], z[order]
        uniq, start = np.unique(lin, return_index=True)
        end = np.r_[start[1:], len(lin)]
        grid.flat[uniq] = np.array([np.nanmax(z[s:e]) for s, e in zip(start, end)])
    else:
        sums = np.zeros_like(grid, dtype=float)
        cnts = np.zeros_like(grid, dtype=int)
        for X, Y, ZZ in zip(xi, yi, z):
            sums[Y, X] += ZZ
            cnts[Y, X] += 1
        with np.errstate(invalid="ignore", divide="ignore"):
            grid = sums / cnts
            grid[cnts == 0] = np.nan
    return grid
# =========================================================================================


# ===================================== MAIN ==============================================
def main():
    datasets = []
    for rid in RUN_IDS:
        x, y, z, meta = load_2d_xy_z(rid)
        x_mV = x * _unit_to_mV(meta["x_unit"])
        y_V  = y * _unit_to_V(meta["y_unit"])
        datasets.append((x_mV, y_V, z, meta))
        print(f"Run {rid} → x:[{np.nanmin(x_mV):.1f},{np.nanmax(x_mV):.1f}] mV, "
              f"y:[{np.nanmin(y_V):.3f},{np.nanmax(y_V):.3f}] V | "
              f"x_name={meta['x_name']}, y_name={meta['y_name']}")

    # y-maximum defined by the chosen run
    y_max_vals = [np.nanmax(y) for (x, y, z, m) in datasets if m["run_id"] == Y_MAX_FROM_RUN]
    if not y_max_vals:
        raise RuntimeError(f"Y_MAX_FROM_RUN={Y_MAX_FROM_RUN} not among RUN_IDS={RUN_IDS}")
    y_max_V = float(y_max_vals[0])

    # Build target grid
    x_edges = np.linspace(X_MIN_MV, X_MAX_MV, NX)
    y_edges = np.linspace(Y_MIN_V,  y_max_V,  NY)

    # Stitch with max
    combined = np.full((NY - 1, NX - 1), np.nan)
    for (x_mV, y_V, z, meta) in datasets:
        g = bin_to_grid(x_mV, y_V, z, x_edges, y_edges, reducer="max")
        if np.isnan(combined).all():
            combined = g
        else:
            both = ~np.isnan(combined) & ~np.isnan(g)
            only_new = np.isnan(combined) & ~np.isnan(g)
            combined[both] = np.maximum(combined[both], g[both])
            combined[only_new] = g[only_new]

    # Plot
    extent = [X_MIN_MV, X_MAX_MV, Y_MIN_V, y_max_V]
    plt.figure(figsize=(7.6, 7.2))
    im = plt.imshow(combined, origin="lower", aspect="auto", extent=extent, interpolation="nearest")
    cbar = plt.colorbar(im, pad=0.02)

    # Labels from metadata
    z_label = datasets[-1][3]["z_label"] or "Signal"
    z_unit  = datasets[-1][3]["z_unit"] or ""
    cbar.set_label(f"{z_label} ({z_unit})" if z_unit else z_label)

    x_label = datasets[-1][3]["x_label"] or "CS gate"
    y_label = datasets[-1][3]["y_label"] or "gate2"
    if "mV" not in x_label: x_label = f"{x_label} (mV)"
    if "V"  not in y_label: y_label = f"{y_label} (V)"

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(f"Stitched runs {RUN_IDS} | x: {X_MIN_MV:.0f}–{X_MAX_MV:.0f} mV, "
              f"y: {Y_MIN_V:.1f}–{y_max_V:.3f} V",
              fontsize=10)
    plt.tight_layout()

    save_path = os.path.join(OUT_DIR, SAVE_NAME)
    plt.savefig(save_path, dpi=300)
    plt.show()
    print(f"✅ Saved: {save_path}")


if __name__ == "__main__":
    main()
# =========================================================================================


