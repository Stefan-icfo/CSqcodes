#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Combine multiple 2D line-sweeps into ONE big plot.

Axes:
- X = CS(inner) in mV (fixed window 760–820 mV)
- Y = main_gate in V (global window −0.3–1.5 V, i.e., "stack the runs vertically")
- Color = G (µS) when available; otherwise the first available dependent.

Cells with no samples remain blank (NaN), so gaps between runs are visible.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ========================= USER SETTINGS =========================
DB_PATH = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v7.db"
RUN_IDS = [809, 810, 811]          # runs to combine

# Global axis windows (exactly as requested)
X_MIN_MV, X_MAX_MV = 760.0, 820.0  # CS(inner) in mV
Y_MIN_V,  Y_MAX_V  = -0.3,  1.5    # main_gate in V

# Grid density (increase for smoother image)
NX, NY = 1001, 901                 # ~0.06 mV horizontal bins; ~2 mV vertical bins
OUT_DIR = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Figures"
OUT_NAME = f"combined_CSx_mainY_runs_{'_'.join(map(str, RUN_IDS))}.png"
# ================================================================

os.makedirs(OUT_DIR, exist_ok=True)
qc.config["core"]["db_location"] = DB_PATH

# --------------------------- helpers ---------------------------
def _flatten(a):
    """Flatten even arrays-of-arrays (snake scans)."""
    arr = np.asarray(a, dtype=object)
    if arr.size == 0:
        return np.array([])
    if arr.dtype != object:
        return np.asarray(arr).ravel()
    return np.concatenate([np.asarray(el).ravel() for el in arr]) if arr.size else np.array([])

def _to_mV_factor(u: str) -> float:
    if not u: return 1.0
    s = u.strip().lower()
    if s in {"mv", "millivolt", "millivolts"}: return 1.0
    if s in {"v", "volt", "volts"}:           return 1e3
    if s in {"uv", "µv", "microvolt", "microvolts"}: return 1e-3
    return 1.0

def _to_V_factor(u: str) -> float:
    if not u: return 1.0
    s = u.strip().lower()
    if s in {"v", "volt", "volts"}:           return 1.0
    if s in {"mv", "millivolt", "millivolts"}: return 1e-3
    if s in {"uv", "µv", "microvolt", "microvolts"}: return 1e-6
    return 1.0

def _label(ds, name, default):
    ps = getattr(ds, "paramspecs", None)
    return getattr(ps.get(name, None), "label", default) if ps else default

def _unit(ds, name, default=""):
    ps = getattr(ds, "paramspecs", None)
    return getattr(ps.get(name, None), "unit", default) if ps else default

def _is_cs(label: str) -> bool:
    s = (label or "").lower()
    return any(t in s for t in ["cs(inner)", "cs (inner)", "cs_inner", "cs", "cs gate", "cs_gate", "csgate", "vcs", "inner"])

def _is_main(label: str) -> bool:
    s = (label or "").lower()
    return any(t in s for t in ["main_gate", "main gate", "main", "outer", "vgo", "gate1", "g1"])

def bin_to_grid(x_mV, y_V, z, x_edges, y_edges):
    """
    Bin scattered (x,y,z) to a regular grid. Cells with no samples stay NaN (blank).
    Aggregation = max per cell (you can change to mean if you prefer).
    """
    xi = np.digitize(x_mV, x_edges) - 1
    yi = np.digitize(y_V,  y_edges) - 1
    nx, ny = len(x_edges) - 1, len(y_edges) - 1

    mask = (xi >= 0) & (xi < nx) & (yi >= 0) & (yi < ny) & np.isfinite(z)
    if not np.any(mask):
        return np.full((ny, nx), np.nan)

    xi, yi, z = xi[mask], yi[mask], z[mask]
    grid = np.full((ny, nx), np.nan)
    lin = yi * nx + xi
    order = np.argsort(lin, kind="mergesort")
    lin, z = lin[order], z[order]
    uniq, start = np.unique(lin, return_index=True)
    end = np.r_[start[1:], len(lin)]
    grid.flat[uniq] = np.array([np.nanmax(z[s:e]) for s, e in zip(start, end)])
    return grid

def extract_csx_mainy(ds):
    """
    From a QCoDeS dataset, pick a 2D block with CS(inner) on X and main_gate on Y.
    Return x (mV), y (V), z (µS or native), and meta dict.
    """
    pdata = ds.get_parameter_data()
    if not pdata:
        raise RuntimeError("empty dataset")

    # prefer G; fall back to I_sens, V_r, Phase, then others
    dep_order = ["G", "I_sens", "V_r", "Phase"] + [k for k in pdata if k not in {"G", "I_sens", "V_r", "Phase"}]

    for dep in dep_order:
        blk = pdata.get(dep, None)
        if blk is None:
            continue
        setpoints = [k for k in blk.keys() if k != dep and not k.endswith("_errors")]
        if len(setpoints) < 2:
            continue

        # We only need one correct pair (first two setpoints usually are)
        a, b = setpoints[0], setpoints[1]
        a_lab = _label(ds, a, a)
        b_lab = _label(ds, b, b)

        # Force CS on X and main_gate on Y (swap if needed)
        if _is_cs(a_lab) and _is_main(b_lab):
            x_name, y_name = a, b
        elif _is_cs(b_lab) and _is_main(a_lab):
            x_name, y_name = b, a
        else:
            # Try to salvage if one looks like CS
            if _is_cs(a_lab):
                x_name, y_name = a, b
            elif _is_cs(b_lab):
                x_name, y_name = b, a
            else:
                continue  # not our pair

        if x_name not in blk or y_name not in blk:
            continue

        x_raw = _flatten(blk[x_name])
        y_raw = _flatten(blk[y_name])
        z_raw = _flatten(blk[dep])
        if x_raw.size == 0 or y_raw.size == 0 or z_raw.size == 0:
            continue

        # Convert x to mV, y to V (IMPORTANT per your request)
        x_mV = x_raw * _to_mV_factor(_unit(ds, x_name, ""))
        y_V  = y_raw * _to_V_factor(_unit(ds, y_name, ""))

        # Convert G to µS if applicable
        z = z_raw.astype(float)
        dep_unit = (_unit(ds, dep, "") or "").lower()
        if dep.lower() == "g" or dep_unit in {"s", "siemens"}:
            z *= 1e6
            cbar_label = "G (µS)"
        else:
            cbar_label = f"{_label(ds, dep, dep)} ({_unit(ds, dep, '')})".strip()

        meta = dict(
            run_id=ds.run_id,
            x_label="CS(inner) (mV)",
            y_label="main_gate (V)",
            cbar_label=cbar_label,
            dep=dep,
            exp_name=getattr(ds, "exp_name", f"run {ds.run_id}")
        )
        return x_mV, y_V, z, meta

    raise RuntimeError("no 2D block with CS on X and main_gate on Y found")

# ----------------------------- main -----------------------------
def main():
    x_edges = np.linspace(X_MIN_MV, X_MAX_MV, NX)
    y_edges = np.linspace(Y_MIN_V,  Y_MAX_V,  NY)

    combined = np.full((NY - 1, NX - 1), np.nan)
    used, skipped = [], []

    for rid in RUN_IDS:
        try:
            ds = qc.load_by_id(rid)
            x_mV, y_V, z, meta = extract_csx_mainy(ds)

            # Keep only samples that fall into the global window
            m = (x_mV >= X_MIN_MV) & (x_mV <= X_MAX_MV) & (y_V >= Y_MIN_V) & (y_V <= Y_MAX_V)
            x_mV, y_V, z = x_mV[m], y_V[m], z[m]

            g = bin_to_grid(x_mV, y_V, z, x_edges, y_edges)

            if np.isnan(combined).all():
                combined = g
            else:
                both = ~np.isnan(combined) & ~np.isnan(g)
                only_new = np.isnan(combined) & ~np.isnan(g)
                combined[both] = np.maximum(combined[both], g[both])
                combined[only_new] = g[only_new]

            used.append(rid)
            print(f"Loaded run {rid}: dep={meta['dep']} → merged.")
        except Exception as e:
            skipped.append((rid, str(e)))
            print(f"Skipping run {rid}: {e}")

    if not used:
        raise RuntimeError("No usable CS×main_gate 2D data found for the selected runs.")

    # Plot one big stitched image
    extent = [X_MIN_MV, X_MAX_MV, Y_MIN_V, Y_MAX_V]
    cmap = plt.cm.magma.copy()     # dark background
    cmap.set_bad(color='black')    # NaN → black

    plt.figure(figsize=(8.8, 8.0))
    im = plt.imshow(combined, origin="lower", aspect="auto",
                    extent=extent, interpolation="nearest", cmap=cmap)
    cbar = plt.colorbar(im, pad=0.02)
    cbar.set_label("G (µS) or chosen signal")

    plt.xlabel("CS(inner) (mV)")
    plt.ylabel("main_gate (V)")
    plt.title(f"Combined runs {used}\nX: {X_MIN_MV:.0f}–{X_MAX_MV:.0f} mV,  Y: {Y_MIN_V:.1f}–{Y_MAX_V:.1f} V")
    plt.tight_layout()

    out = os.path.join(OUT_DIR, OUT_NAME)
    plt.savefig(out, dpi=300)
    plt.show()
    print("Saved:", out)

    if skipped:
        print("\nSkipped runs:")
        for rid, msg in skipped:
            print(f"  - {rid}: {msg}")

if __name__ == "__main__":
    main()





