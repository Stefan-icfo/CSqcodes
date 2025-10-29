#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ===================== USER CONFIG =====================
run_data = 288
run_background = 282
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v19_211025.db"
out_dir  = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Figures"
os.makedirs(out_dir, exist_ok=True)

# Subtraction mode:
DO_VOLTAGE_DOMAIN = True   # True: V_sig = sqrt(P_data) - sqrt(P_bg);  P_sig = V_sig**2
USE_INDEX_AXIS    = False  # MUST be False for physical integration

# Averaging (plot/also for averaged-area):
AVG_BLOCK = 4  # non-overlapping block size

# ---- Peak integration window (MHz) ----
PEAK_START_MHZ = 152.9671
PEAK_END_MHZ   = 152.97092

# ---- Baseline removal inside the window ('none' | 'constant' | 'linear') ----
BASELINE_MODE = 'none'
BASELINE_EDGE_FRACTION = 0.1  # fraction of points at the edges used for baseline estimation

# ---- Fill toggles ----
# NOTE: The RED fill section is commented out below. To enable it later, remove the leading '#'.
FILL_AVG_GREEN = True
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
    N = len(y_window)
    if N < 4 or mode not in ('constant', 'linear'):
        return np.zeros_like(y_window)

    k = max(1, int(np.floor(edge_frac * N)))
    left_idx  = np.arange(k)
    right_idx = np.arange(N - k, N)
    idx = np.concatenate([left_idx, right_idx])
    x = np.arange(N, dtype=float)

    if mode == 'constant':
        b = np.mean(y_window[idx])
        return np.full_like(y_window, b)

    # linear
    A = np.vstack([x[idx], np.ones_like(idx, dtype=float)]).T
    a, b = np.linalg.lstsq(A, y_window[idx], rcond=None)[0]  # y ~ a*x + b
    return a * x + b

def integrate_trapz(x_Hz, y_W_per_Hz):
    """Trapezoidal integration ∫ y(f) df over Hz -> W."""
    order = np.argsort(x_Hz)
    return np.trapz(y_W_per_Hz[order], x_Hz[order])

# --------------- checks ---------------
if USE_INDEX_AXIS:
    raise RuntimeError("Physical integration requires frequency axis. Set USE_INDEX_AXIS = False.")
if not (PEAK_START_MHZ < PEAK_END_MHZ):
    raise ValueError("PEAK_START_MHZ must be < PEAK_END_MHZ.")

# --------------- load & subtract ---------------
x_d, y_d, dep_d, fkey_d = load_1d_in_order(run_data)
x_b, y_b, dep_b, fkey_b = load_1d_in_order(run_background)

n = min(len(y_d), len(y_b))
if n == 0:
    raise RuntimeError("No points to subtract (one of the traces is empty).")

x_MHz = x_d[:n]  # frequency axis (MHz)
if fkey_d is None:
    raise RuntimeError("No frequency key found: cannot integrate.")

if DO_VOLTAGE_DOMAIN:
    V_d = np.sqrt(np.clip(y_d[:n], 0.0, None))
    V_b = np.sqrt(np.clip(y_b[:n], 0.0, None))
    V_sig = V_d - V_b
    y_sig = V_sig**2             # W/Hz
    y_label = "PSD (pW/Hz)"
else:
    y_sig = y_d[:n] - y_b[:n]    # W/Hz
    y_label = "PSD (W/Hz)  [power subtraction]"

print(f"[subtract] using first {n} points (index-by-index).")
print(f"  data range: [{np.nanmin(y_d[:n]):.3e}, {np.nanmax(y_d[:n]):.3e}]")
print(f"  back range: [{np.nanmin(y_b[:n]):.3e}, {np.nanmax(y_b[:n]):.3e}]")
print(f"  result rng: [{np.nanmin(y_sig):.3e}, {np.nanmax(y_sig):.3e}]")

# For display (pW/Hz)
y_sig_plot = y_sig / 1e-12

# --------------- averaging (used for both plot and averaged-area) ---------------
x_avg, y_avg = block_average(x_MHz, y_sig_plot, AVG_BLOCK)

# --------------- RAW integration window ---------------
f0, f1 = float(PEAK_START_MHZ), float(PEAK_END_MHZ)
mask_raw = (x_MHz >= f0) & (x_MHz <= f1)
if not np.any(mask_raw):
    raise RuntimeError("No points inside the selected [PEAK_START_MHZ, PEAK_END_MHZ] window (raw).")

x_win_raw_MHz = x_MHz[mask_raw]
y_win_raw_W_per_Hz = y_sig[mask_raw]  # SI units for integration

# optional baseline inside window (raw)
if BASELINE_MODE.lower() in ('constant', 'linear'):
    y_win_raw_pW = y_win_raw_W_per_Hz / 1e-12
    y_base_raw = estimate_baseline(y_win_raw_pW, mode=BASELINE_MODE.lower(), edge_frac=BASELINE_EDGE_FRACTION)
    y_win_raw_pW_corr = y_win_raw_pW - y_base_raw
    y_win_raw_W_per_Hz = y_win_raw_pW_corr * 1e-12  # back to SI

# integrate RAW
x_win_raw_Hz = x_win_raw_MHz * 1e6
area_raw_W = integrate_trapz(x_win_raw_Hz, y_win_raw_W_per_Hz)
area_raw_pW = area_raw_W * 1e12

# --------------- AVERAGED integration window ---------------
area_avg_W = None
area_avg_pW = None
x_win_avg_MHz = None
y_win_avg_W_per_Hz = None
if x_avg.size and y_avg.size:
    mask_avg = (x_avg >= f0) & (x_avg <= f1)
    if np.any(mask_avg):
        x_win_avg_MHz = x_avg[mask_avg]
        y_win_avg_pW_per_Hz = y_avg[mask_avg]  # already in pW/Hz
        # optional baseline on averaged
        if BASELINE_MODE.lower() in ('constant', 'linear'):
            y_base_avg = estimate_baseline(y_win_avg_pW_per_Hz, mode=BASELINE_MODE.lower(), edge_frac=BASELINE_EDGE_FRACTION)
            y_win_avg_pW_per_Hz = y_win_avg_pW_per_Hz - y_base_avg
        # integrate averaged
        x_win_avg_Hz = x_win_avg_MHz * 1e6
        y_win_avg_W_per_Hz = y_win_avg_pW_per_Hz * 1e-12
        area_avg_W = integrate_trapz(x_win_avg_Hz, y_win_avg_W_per_Hz)
        area_avg_pW = area_avg_W * 1e12
    else:
        print("Warning: no averaged points fall inside the integration window.")
else:
    print("Warning: averaged arrays empty (check AVG_BLOCK).")

print("\n=== Peak Area (noise-subtracted) ===")
print(f"Window: [{f0:.6f}, {f1:.6f}] MHz   ({f1 - f0:.6f} MHz wide) | Baseline: {BASELINE_MODE}")
print(f"RAW      area = {area_raw_W:.6e} W   = {area_raw_pW:.6e} pW")
if area_avg_W is not None:
    print(f"AVERAGED area = {area_avg_W:.6e} W   = {area_avg_pW:.6e} pW")

# --------------- plot ---------------
fig, ax = plt.subplots(figsize=(10, 6))

# averaged curve for visibility (points + line)
if x_avg.size and y_avg.size:
    ax.plot(x_avg, y_avg, '-o', ms=3.5, lw=1.8, label=f"{AVG_BLOCK}-point average")

# --- RED FILL (currently disabled) ---
# The block below would fill in RED the RAW integrated area.
# It is commented out on purpose so the red area is NOT shown now.
# If you want to show it in the future, remove the leading '#' on each line.

# x_fill_raw = x_win_raw_MHz.copy()
# y_fill_raw_pW = (y_win_raw_W_per_Hz / 1e-12).copy()
# order = np.argsort(x_fill_raw)
# x_fill_raw = x_fill_raw[order]
# y_fill_raw_pW = y_fill_raw_pW[order]
# y_fill_raw_pW = np.clip(y_fill_raw_pW, 0.0, None)  # visual clipping only
# ax.fill_between(x_fill_raw, 0.0, y_fill_raw_pW, color='crimson', alpha=0.35,
#                 label='Integrated area (RAW, red)')

# --- GREEN FILL: averaged integrated area (shown) ---
if FILL_AVG_GREEN and (x_win_avg_MHz is not None) and (y_win_avg_W_per_Hz is not None):
    x_fill_avg = x_win_avg_MHz.copy()
    y_fill_avg_pW = (y_win_avg_W_per_Hz / 1e-12).copy()
    order = np.argsort(x_fill_avg)
    x_fill_avg = x_fill_avg[order]
    y_fill_avg_pW = y_fill_avg_pW[order]
    y_fill_avg_pW = np.clip(y_fill_avg_pW, 0.0, None)
    ax.fill_between(x_fill_avg, 0.0, y_fill_avg_pW, color='green', alpha=0.30,
                    label='Integrated area (AVG, green)')

ax.set_xlabel("Frequency (MHz)", fontsize=18, fontname="Calibri")
ax.set_ylabel(y_label, fontsize=18, fontname="Calibri")
ax.set_title("Noise-subtracted PSD: integrated areas (AVG shown in green)", fontsize=18, fontname="Calibri")

ax.tick_params(axis='both', which='major', labelsize=16)
for lab in (ax.get_xticklabels() + ax.get_yticklabels()):
    lab.set_fontname("Calibri")

# Annotation box with both area values
annot = f"RAW area = {area_raw_pW:.3e} pW"
if area_avg_pW is not None:
    annot += f"\nAVG area = {area_avg_pW:.3e} pW"
annot += f"\n[{PEAK_START_MHZ:.6f}, {PEAK_END_MHZ:.6f}] MHz, baseline: {BASELINE_MODE}"
ax.text(0.02, 0.98, annot, transform=ax.transAxes, ha='left', va='top')

ax.legend(frameon=False)
out_png = os.path.join(out_dir, f"psd_area_avg_green_raw_red_off_run{run_data}_minus_{run_background}_avg{AVG_BLOCK}.png")
plt.savefig(out_png, dpi=220, bbox_inches='tight')
print("Saved figure to:", out_png)

try:
    plt.show()
except Exception as e:
    print("plt.show() failed:", e)


