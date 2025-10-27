#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ===================== USER CONFIG =====================
run_data =  4088
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v11.db"
out_dir  = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Figures"
os.makedirs(out_dir, exist_ok=True)

# Averaging (only averaged curve will be plotted)
AVG_BLOCK = 9     # set 1 to disable averaging (plot will be identical to raw)

# ---- Peak integration window (MHz) ----
PEAK_START_MHZ = 153.65823
PEAK_END_MHZ   = 153.66190

# ---- Baseline removal inside the window ('none' | 'constant' | 'linear') ----
BASELINE_MODE = 'none'
BASELINE_EDGE_FRACTION = 0.10

# ---- Plot style ----
BORDEAUX_PASTEL = "#B04A5A"   # averaged curve color (pastel bordeaux)
# =======================================================

qc.config["core"]["db_location"] = db_path

# --------------- helpers ---------------
CANDIDATE_SIGNAL_KEYS = [
    "avg_avg_psd_nodrive", "avg_psd_nodrive",
    "avg_avg_psd", "avg_psd",
    "V_fft_avg_avg", "PSD", "psd"
]
FREQ_KEY_CANDIDATES = ["freq_param", "frequency", "f", "freq"]

def load_1d_in_order(run_id):
    """Load a 1D trace in ACQUISITION ORDER (no sorting). Returns (x_MHz, y, dep_key, freq_key)."""
    ds = qc.load_by_id(run_id)
    pd = ds.get_parameter_data()

    dep = next((k for k in CANDIDATE_SIGNAL_KEYS if k in pd), next(iter(pd.keys())))
    blk = pd[dep]

    freq_key = None
    for k in blk.keys():
        if k == dep:
            continue
        if any(t in k.lower() for t in FREQ_KEY_CANDIDATES):
            freq_key = k
            break
    if freq_key is None:
        freq_key = next((k for k in blk.keys() if k != dep), None)

    x = np.asarray(blk[freq_key]).ravel() / 1e6 if freq_key else np.arange(len(np.asarray(blk[dep]).ravel()), dtype=float)
    y = np.asarray(blk[dep]).ravel()

    m = np.isfinite(y)
    if freq_key is not None:
        m &= np.isfinite(x)
        x = x[m]
    else:
        x = np.arange(np.count_nonzero(m), dtype=float)
    y = y[m]

    print(f"[run {run_id}] dep='{dep}', f_key='{freq_key}', N={len(y)}")
    return x, y, dep, freq_key

def block_average(x, y, n):
    """Non-overlapping block average of size n; tail truncated."""
    if n <= 1:
        return x.copy(), y.copy()
    m = (len(y) // n) * n
    if m == 0:
        return np.array([]), np.array([])
    xb = x[:m].reshape(-1, n).mean(axis=1)
    yb = y[:m].reshape(-1, n).mean(axis=1)
    return xb, yb

def estimate_baseline(y_window, mode='none', edge_frac=0.1):
    """
    Baseline inside the window.
    - 'none'     : zeros
    - 'constant' : mean of edges (left+right)
    - 'linear'   : least-squares line from edges
    """
    mode = (mode or 'none').lower()
    N = len(y_window)
    if mode not in ('constant', 'linear') or N < 4:
        return np.zeros_like(y_window)

    k = max(1, int(np.floor(edge_frac * N)))
    left_idx  = np.arange(k)
    right_idx = np.arange(N - k, N)
    idx = np.concatenate([left_idx, right_idx])
    x = np.arange(N, dtype=float)

    if mode == 'constant':
        b = np.mean(y_window[idx])
        return np.full_like(y_window, b)

    A = np.vstack([x[idx], np.ones_like(idx, dtype=float)]).T
    a, b = np.linalg.lstsq(A, y_window[idx], rcond=None)[0]  # y ~ a*x + b
    return a * x + b

def integrate_trapz(x_Hz, y_W_per_Hz):
    """Trapezoidal integration ∫ y(f) df over Hz -> W."""
    order = np.argsort(x_Hz)
    return np.trapz(y_W_per_Hz[order], x_Hz[order])

# --------------- checks ---------------
if not (PEAK_START_MHZ < PEAK_END_MHZ):
    raise ValueError("PEAK_START_MHZ must be < PEAK_END_MHZ.")

# --------------- load (NO background subtraction) ---------------
x_MHz, y_psd, dep_key, fkey = load_1d_in_order(run_data)
if fkey is None:
    raise RuntimeError("No frequency key found: cannot integrate physically.")

# y_psd assumed in SI units (W/Hz). For display use pW/Hz.
y_psd_pW_per_Hz = y_psd / 1e-12

# --------------- averaging (ONLY this will be plotted) ---------------
x_avg, y_avg_pW_per_Hz = block_average(x_MHz, y_psd_pW_per_Hz, AVG_BLOCK)

# --------------- Integration window (optional, not shown on plot) ---------------
f0, f1 = float(PEAK_START_MHZ), float(PEAK_END_MHZ)
mask_raw = (x_MHz >= f0) & (x_MHz <= f1)
if not np.any(mask_raw):
    raise RuntimeError("No points inside the selected [PEAK_START_MHZ, PEAK_END_MHZ] window.")

x_win_raw_MHz = x_MHz[mask_raw]
y_win_raw_pW_per_Hz = y_psd_pW_per_Hz[mask_raw]

if BASELINE_MODE.lower() in ('constant', 'linear'):
    y_base_raw = estimate_baseline(y_win_raw_pW_per_Hz, mode=BASELINE_MODE, edge_frac=BASELINE_EDGE_FRACTION)
    y_win_raw_pW_per_Hz = y_win_raw_pW_per_Hz - y_base_raw

# integrate (values printed only to console)
x_win_raw_Hz = x_win_raw_MHz * 1e6
area_raw_W = integrate_trapz(x_win_raw_Hz, y_win_raw_pW_per_Hz * 1e-12)
area_raw_pW = area_raw_W * 1e12
print("\n=== Peak Area (no background subtraction; NOT shown on plot) ===")
print(f"Window: [{f0:.6f}, {f1:.6f}] MHz   ({f1 - f0:.6f} MHz wide) | Baseline: {BASELINE_MODE}")
print(f"AREA = {area_raw_W:.6e} W   = {area_raw_pW:.6e} pW")

# --------------- plot (ONLY averaged curve, bordeaux pastel) ---------------
fig, ax = plt.subplots(figsize=(10, 6))

if x_avg.size and y_avg_pW_per_Hz.size:
    ax.plot(
        x_avg, y_avg_pW_per_Hz,
        '-o', ms=3.5, lw=1.8,
        label=f"{AVG_BLOCK}-point average",
        color=BORDEAUX_PASTEL
    )
else:
    # fallback: plot raw (still in bordeaux) only if averaging is empty
    ax.plot(x_MHz, y_psd_pW_per_Hz, lw=1.2, color=BORDEAUX_PASTEL, label="PSD")

ax.set_xlabel("Frequency +153.6 (MHz)", fontsize=18, fontname="Calibri")
ax.set_ylabel("PSD (pW/Hz)", fontsize=18, fontname="Calibri")
ax.set_title("PSD (no background): averaged trace", fontsize=18, fontname="Calibri")

ax.tick_params(axis='both', which='major', labelsize=16)
for lab in (ax.get_xticklabels() + ax.get_yticklabels()):
    try:
        lab.set_fontname("Calibri")
    except Exception:
        pass

ax.legend(frameon=False)
out_png = os.path.join(out_dir, f"psd_avg_bordeaux_NO_BG_run{run_data}_avg{AVG_BLOCK}.png")
plt.savefig(out_png, dpi=220, bbox_inches='tight')
print("Saved figure to:", out_png)

try:
    plt.show()
except Exception as e:
    print("plt.show() failed:", e)




