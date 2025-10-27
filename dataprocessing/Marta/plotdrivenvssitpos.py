#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Driven frequency-sweep → 2D map (no background subtraction, no unit conversion).
X = frequency (MHz)
Y = CS gate parsed from exp_name ('gcs=...')
Color = driven intensity (auto-detected dependent; raw units)

- Auto-detects frequency key (e.g., 'zurich_oscs0_freq', 'freq_param', ...).
- Auto-detects a sensible dependent channel (V_r → V_fft_avg_avg → Amplitude → ... → avg_avg_psd_nodrive).
- Interpolates all traces to a common frequency grid (intersection).
- Sorts rows by gate voltage. Handles a single gate gracefully.

Set RUN_IDS and DB path below.
"""

import re
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ===================== USER SETTINGS =====================
qc.config["core"]["db_location"] = (
    r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v11.db"
)

RUN_IDS = list(range(2252, 2307, 2))   # your driven runs (no background run needed)

# (Optional) force keys. Leave None to auto-detect.
FORCE_DEP_KEY  = None      # e.g. "V_r" or "V_fft_avg_avg"
FORCE_FREQ_KEY = None      # e.g. "zurich_oscs0_freq"

# Plot presentation
TITLE = "Driven map: frequency vs CS gate"
NORMALIZE_PER_TRACE = False     # True → divide each row by its max
USE_LOG_COLOR = False           # True → log10 color scale
CMAP = "viridis"
# =========================================================

# Preferred dependents to try (in order) when FORCE_DEP_KEY is None:
PREF_DEP_KEYS = [
    "V_r", "V_fft_avg_avg", "Amplitude", "R", "I_rf", "I", "Phase",
    "avg_avg_psd_nodrive", "avg_psd_nodrive", "avg_avg_psd"
]

def extract_gcs_voltage(name: str):
    """Parse 'gcs=...' from exp_name; returns float or None."""
    m = re.search(r"gcs=([+-]?\d+(?:\.\d+)?)", str(name))
    return float(m.group(1)) if m else None

def to_mV_if_volts(arr):
    """Convert gate axis to mV if values look like volts (<~5 abs)."""
    a = np.asarray(arr, dtype=float)
    if a.size == 0:
        return a, "arb."
    return (a * 1e3, "mV") if np.nanmax(np.abs(a)) < 5.0 else (a, "V")

def pick_frequency_key(block, dep):
    """
    Choose a frequency-like setpoint key in this parameter block.
    Priority: FORCE_FREQ_KEY → any setpoint containing 'freq' → first other setpoint.
    """
    if FORCE_FREQ_KEY and FORCE_FREQ_KEY in block:
        return FORCE_FREQ_KEY
    candidates = [k for k in block.keys() if k != dep and not k.endswith("_errors")]
    for k in candidates:
        if "freq" in k.lower():
            return k
    if candidates:
        return candidates[0]
    raise KeyError("No suitable frequency setpoint key found in this block.")

def pick_dependent_key(param_data):
    """
    Choose a dependent present in the dataset.
    Priority: FORCE_DEP_KEY → in PREF_DEP_KEYS order → first available key.
    """
    keys = list(param_data.keys())
    if FORCE_DEP_KEY:
        if FORCE_DEP_KEY in param_data:
            return FORCE_DEP_KEY
        raise KeyError(f"Forced dependent '{FORCE_DEP_KEY}' not found. Available: {keys}")
    for k in PREF_DEP_KEYS:
        if k in param_data:
            return k
    if not keys:
        raise KeyError("Dataset has no parameters.")
    return keys[0]

def load_one_run(run_id):
    """
    Load one driven sweep from QCoDeS:
      returns (f_MHz_sorted, y_sorted, dep_label, dep_unit, ds)
    """
    ds = qc.load_by_id(run_id)
    pd = ds.get_parameter_data()
    if not pd:
        raise RuntimeError(f"Run {run_id}: empty dataset")

    dep = pick_dependent_key(pd)
    block = pd[dep]
    fkey = pick_frequency_key(block, dep)

    y = np.asarray(block[dep]).ravel().astype(float)
    f = np.asarray(block[fkey]).ravel().astype(float) / 1e6  # Hz → MHz

    m = np.isfinite(f) & np.isfinite(y)
    if not np.any(m):
        raise ValueError(f"Run {run_id}: no finite points")
    f, y = f[m], y[m]

    # sort by frequency
    idx = np.argsort(f)
    f, y = f[idx], y[idx]

    # labels/units if present
    ps = getattr(ds, "paramspecs", None)
    dep_label = dep
    dep_unit = ""
    if ps and ps.get(dep):
        dep_label = getattr(ps[dep], "label", dep_label) or dep_label
        dep_unit  = getattr(ps[dep], "unit", dep_unit) or dep_unit

    return f, y, dep_label, dep_unit, ds

def safe_extent(x, y, xname="x", yname="y"):
    x = np.asarray(x); y = np.asarray(y)
    if x.size == 0: raise ValueError(f"{xname} axis is empty")
    if y.size == 0: raise ValueError(f"{yname} axis is empty")
    # Handle single-row nicely: pad a little so imshow has non-zero height
    ymin, ymax = np.nanmin(y), np.nanmax(y)
    if ymax == ymin:
        pad = max(1e-3 * max(1.0, abs(ymin)), 1e-6)
        ymin -= pad; ymax += pad
    return [np.nanmin(x), np.nanmax(x), ymin, ymax]

# --------- load all runs ---------
xs, ys, gates = [], [], []
kept, skipped = [], []
dep_label_used, dep_unit_used = None, None

for rid in RUN_IDS:
    try:
        f, val, dep_label, dep_unit, ds = load_one_run(rid)
        vg = extract_gcs_voltage(getattr(ds, "exp_name", ""))
        if vg is None:
            skipped.append((rid, "exp_name has no 'gcs='"))
            continue
        xs.append(f); ys.append(val); gates.append(vg); kept.append(rid)
        if dep_label_used is None:
            dep_label_used, dep_unit_used = dep_label, dep_unit
    except Exception as e:
        skipped.append((rid, str(e)))

if not ys:
    raise RuntimeError(
        "No usable runs.\n" + "\n".join([f"- {rid}: {msg}" for rid, msg in skipped])
    )

# --------- common frequency grid (intersection) ---------
fmins = [np.nanmin(f) for f in xs]
fmaxs = [np.nanmax(f) for f in xs]
fmin = float(np.max(fmins))
fmax = float(np.min(fmaxs))
if not np.isfinite(fmin) or not np.isfinite(fmax) or fmax <= fmin:
    raise RuntimeError("Empty frequency intersection across runs.")

n_common = int(np.min([len(f) for f in xs]))  # safe number of columns
freqs = np.linspace(fmin, fmax, n_common)

# Interpolate to common grid
mat = np.vstack([np.interp(freqs, f, y) for f, y in zip(xs, ys)])

# Optional per-trace normalization
if NORMALIZE_PER_TRACE:
    s = np.nanmax(np.abs(mat), axis=1, keepdims=True)
    s[s == 0] = 1.0
    mat = mat / s

# Gate axis and sorting
gates = np.asarray(gates, dtype=float)
gates_plot, gate_unit = to_mV_if_volts(gates)
order = np.argsort(gates_plot)
mat = mat[order]; gates_plot = gates_plot[order]

# --------- plot ---------
extent = safe_extent(freqs, gates_plot, "freqs (MHz)", f"V_gcs ({gate_unit})")

plt.figure(figsize=(10, 6))
if USE_LOG_COLOR:
    eps = np.nanmax(mat) * 1e-9 if np.nanmax(mat) > 0 else 1e-12
    data = np.log10(np.clip(mat, eps, None))
    im = plt.imshow(data, aspect="auto", origin="lower", extent=extent,
                    cmap=CMAP, interpolation="nearest")
    cbar_label = f"log10({dep_label_used} [{dep_unit_used}])" if dep_unit_used else f"log10({dep_label_used})"
else:
    im = plt.imshow(mat, aspect="auto", origin="lower", extent=extent,
                    cmap=CMAP, interpolation="nearest")
    cbar_label = f"{dep_label_used} [{dep_unit_used}]" if dep_unit_used else f"{dep_label_used}"

cbar = plt.colorbar(im, pad=0.02)
cbar.set_label(cbar_label)

plt.xlabel("Frequency (MHz)")
plt.ylabel(f"V_gcs ({gate_unit})")
title_runs = f"runs {kept[0]}..{kept[-1]}" if len(kept) > 1 else f"run {kept[0]}"
plt.title(f"{TITLE} — {title_runs} (kept {len(kept)}, skipped {len(skipped)})")
plt.tight_layout()
plt.show()

if skipped:
    print("\nSkipped runs:")
    for rid, msg in skipped:
        print(f"  - {rid}: {msg}")




