#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Single-run 2D plot with auto-ranging:
- X = CS(inner) in mV (auto min/max from data)
- Y = main_gate in V (auto min/max from data)
- Color = G (uS) if present, else first available dependent.
No manual X/Y limits required.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ======= USER SETTINGS =======
DB_PATH = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v11.db"
RUN_ID  = 1914                      # <-- your single run id
NX, NY  = 1001, 901                # grid density (bigger = smoother)
PADDING_FRAC = 0.02                # 2% padding around data range
OUT_DIR = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Figures"
# ==============================

os.makedirs(OUT_DIR, exist_ok=True)
qc.config["core"]["db_location"] = DB_PATH

# ---------- utilities ----------
def _flatten(a):
    arr = np.asarray(a, dtype=object)
    if arr.size == 0: return np.array([])
    if arr.dtype != object: return np.asarray(arr).ravel()
    return np.concatenate([np.asarray(el).ravel() for el in arr]) if arr.size else np.array([])

def _unit(ds, name, default=""):
    ps = getattr(ds, "paramspecs", None)
    return getattr(ps.get(name, None), "unit", default) if ps else default

def _label(ds, name, default):
    ps = getattr(ds, "paramspecs", None)
    return getattr(ps.get(name, None), "label", default) if ps else default

def _is_cs(label: str) -> bool:
    s = (label or "").lower()
    return any(t in s for t in ["cs(inner)", "cs (inner)", "cs_inner", "cs", "cs gate", "csgate", "vcs", "inner"])

def _is_main(label: str) -> bool:
    s = (label or "").lower()
    return any(t in s for t in ["main_gate", "main gate", "main", "outer", "vgo", "gate1", "g1"])

def _to_mV_factor(u: str) -> float:
    s = (u or "").strip().lower()
    if s in {"mv","millivolt","millivolts",""}: return 1.0
    if s in {"v","volt","volts"}:               return 1e3
    if s in {"uv","µv","microvolt","microvolts"}: return 1e-3
    return 1.0

def _to_V_factor(u: str) -> float:
    s = (u or "").strip().lower()
    if s in {"v","volt","volts",""}:            return 1.0
    if s in {"mv","millivolt","millivolts"}:     return 1e-3
    if s in {"uv","µv","microvolt","microvolts"}: return 1e-6
    return 1.0

def bin_to_grid(x_mV, y_V, z, x_edges, y_edges):
    """Bin scattered (x,y,z) onto a regular grid with max aggregator; NaNs remain blank."""
    xi = np.digitize(x_mV, x_edges) - 1
    yi = np.digitize(y_V,  y_edges) - 1
    nx, ny = len(x_edges) - 1, len(y_edges) - 1
    mask = (xi >= 0) & (xi < nx) & (yi >= 0) & (yi < ny) & np.isfinite(z)
    grid = np.full((ny, nx), np.nan)
    if not np.any(mask): return grid
    xi, yi, z = xi[mask], yi[mask], z[mask]
    lin = yi * nx + xi
    order = np.argsort(lin, kind="mergesort")
    lin, z = lin[order], z[order]
    uniq, start = np.unique(lin, return_index=True)
    end = np.r_[start[1:], len(lin)]
    grid.flat[uniq] = np.array([np.nanmax(z[s:e]) for s, e in zip(start, end)])
    return grid

def extract_one(ds):
    """Pick a 2D block with CS(inner) on X, main_gate on Y. Return x(mV), y(V), z, labels."""
    pdata = ds.get_parameter_data()
    if not pdata: raise RuntimeError("Empty dataset")

    prefer = ["G", "I_sens", "V_r", "Phase"] + [k for k in pdata if k not in {"G","I_sens","V_r","Phase"}]
    for dep in prefer:
        blk = pdata.get(dep)
        if blk is None: continue
        setpoints = [k for k in blk.keys() if k != dep and not k.endswith("_errors")]
        if len(setpoints) < 2: continue

        a, b = setpoints[0], setpoints[1]
        La, Lb = _label(ds, a, a), _label(ds, b, b)

        if _is_cs(La) and _is_main(Lb): x_name, y_name = a, b
        elif _is_cs(Lb) and _is_main(La): x_name, y_name = b, a
        else: continue

        x_raw = _flatten(blk[x_name]); y_raw = _flatten(blk[y_name]); z_raw = _flatten(blk[dep])
        if x_raw.size == 0 or y_raw.size == 0 or z_raw.size == 0: continue

        x_mV = x_raw * _to_mV_factor(_unit(ds, x_name, ""))
        y_V  = y_raw * _to_V_factor(_unit(ds, y_name, ""))

        z = z_raw.astype(float)
        dep_unit = (_unit(ds, dep, "") or "").lower()
        if dep.lower()=="g" or dep_unit in {"s","siemens"}:
            z *= 1e6
            cbar = "G (uS)"
        else:
            label = _label(ds, dep, dep)
            unit  = _unit(ds, dep, "")
            cbar  = f"{label} ({unit})".strip()

        return x_mV, y_V, z, cbar

    raise RuntimeError("No 2D CS(inner) x main_gate block found")

# ---------- main ----------
def main():
    ds = qc.load_by_id(RUN_ID)
    x_mV, y_V, z, cbar_label = extract_one(ds)

    # auto-range with small padding
    xmin, xmax = float(np.nanmin(x_mV)), float(np.nanmax(x_mV))
    ymin, ymax = float(np.nanmin(y_V)),  float(np.nanmax(y_V))
    dx = (xmax - xmin) or 1.0
    dy = (ymax - ymin) or 1.0
    xmin -= PADDING_FRAC * dx; xmax += PADDING_FRAC * dx
    ymin -= PADDING_FRAC * dy; ymax += PADDING_FRAC * dy

    x_edges = np.linspace(xmin, xmax, NX)
    y_edges = np.linspace(ymin, ymax, NY)

    grid = bin_to_grid(x_mV, y_V, z, x_edges, y_edges)

    extent = [xmin, xmax, ymin, ymax]
    cmap = plt.cm.magma.copy(); cmap.set_bad('black')

    plt.figure(figsize=(8.2, 7.4))
    im = plt.imshow(grid, origin="lower", aspect="auto", extent=extent,
                    cmap=cmap, interpolation="nearest")
    cbar = plt.colorbar(im, pad=0.02); cbar.set_label(cbar_label)
    plt.xlabel("CS(inner) (mV)")
    plt.ylabel("main_gate (V)")
    plt.title(f"Run {RUN_ID} (auto range)")
    plt.tight_layout()

    out = os.path.join(OUT_DIR, f"single_run_{RUN_ID}_auto.png")
    plt.savefig(out, dpi=300); plt.show()
    print("Saved:", out)

if __name__ == "__main__":
    main()
