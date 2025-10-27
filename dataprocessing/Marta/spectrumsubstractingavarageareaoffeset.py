#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ===================== USER CONFIG =====================
# ===================== USER CONFIG =====================
run_data = 92
run_background = 83
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v13.db"
out_dir  = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Figures"
os.makedirs(out_dir, exist_ok=True)
# Subtraction mode:
DO_VOLTAGE_DOMAIN = True   # True: V_sig = sqrt(P_data) - sqrt(P_bg);  P_sig = V_sig**2
USE_INDEX_AXIS    = False  # MUST be False for physical integration

# Averaging (plot/also for averaged-area):
AVG_BLOCK = 5         # non-overlapping block size

# ---- Peak integration window (MHz) ----
PEAK_START_MHZ = 153.336
PEAK_END_MHZ   = 153.3388

# ---- Baseline options inside the window:
# 'none' | 'constant' | 'linear' | 'manual'
BASELINE_MODE = 'manual'
BASELINE_EDGE_FRACTION = 0.10    # edges used for 'constant'/'linear'
MANUAL_BASELINE_PW_PER_HZ =0# <-- set your horizontal line here (in pW/Hz)

# ---- Fill toggles ----
# (RED RAW area re
# mains commented out below; remove '#' to show it.)
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

def estimate_baseline(y_window_pW, mode='none', edge_frac=0.1, manual_pW=0):
    """
    Return the baseline array (in pW/Hz) of the same length as y_window_pW.
    - 'none'     : zeros
    - 'manual'   : flat line at 'manual_pW'
    - 'constant' : mean of left+right edges
    - 'linear'   : least-squares line from edges
    """
    N = len(y_window_pW)
    if mode == 'none' or N < 4:
        return np.zeros_like(y_window_pW)
    if mode == 'manual':
        return np.full_like(y_window_pW, float(manual_pW))
    if mode not in ('constant', 'linear'):
        return np.zeros_like(y_window_pW)

    k = max(1, int(np.floor(edge_frac * N)))
    left_idx  = np.arange(k)
    right_idx = np.arange(N - k, N)
    idx = np.concatenate([left_idx, right_idx])
    x = np.arange(N, dtype=float)

    if mode == 'constant':
        b = np.mean(y_window_pW[idx])
        return np.full_like(y_window_pW, b)

    # linear
    A = np.vstack([x[idx], np.ones_like(idx, dtype=float)]).T
    a, b = np.linalg.lstsq(A, y_window_pW[idx], rcond=None)[0]  # y ~ a*x + b
    return a * x + b

def integrate_trapz(x_Hz, y_W_per_Hz):
    """Trapezoidal integration ∫ y(f) df over Hz -> W."""
    order = np.argsort(x_Hz)
    return np.trapz(y_W_per_Hz[order], x_Hz[order])

# --------------- checks ---------------
if USE_INDEX_AXIS:
    raise RuntimeError("Physical integration requires frequency axis. Set USE_INDEX_AXIS = False.")

# --------------- load & subtract ---------------
x_d, y_d, dep_d, fkey_d = load_1d_in_order(run_data)
x_b, y_b, dep_b, fkey_b = load_1d_in_order(run_background)

n = min(len(y_d), len(y_b))
if n == 0:
    raise RuntimeError("No points to subtract (one trace empty).")

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

print(f"[subtract] using first {n} points.")
print(f"  data range: [{np.nanmin(y_d[:n]):.3e}, {np.nanmax(y_d[:n]):.3e}]")
print(f"  back range: [{np.nanmin(y_b[:n]):.3e}, {np.nanmax(y_b[:n]):.3e}]")
print(f"  result rng: [{np.nanmin(y_sig):.3e}, {np.nanmax(y_sig):.3e}]")

# For display (pW/Hz)
y_sig_plot = y_sig / 1e-12

# --------------- averaging (plot + averaged-area) ---------------
x_avg, y_avg = block_average(x_MHz, y_sig_plot, AVG_BLOCK)

# --------------- integration window ---------------
f0, f1 = float(PEAK_START_MHZ), float(PEAK_END_MHZ)
mask_raw = (x_MHz >= f0) & (x_MHz <= f1)
if not np.any(mask_raw):
    raise RuntimeError("No points inside the selected window.")

x_win_raw_MHz = x_MHz[mask_raw]
y_win_raw_W_per_Hz = y_sig[mask_raw]        # SI

# ----- baseline (in pW/Hz), then integrate only above baseline -----
y_win_raw_pW = y_win_raw_W_per_Hz / 1e-12
baseline_raw_pW = estimate_baseline(
    y_win_raw_pW, mode=BASELINE_MODE.lower(),
    edge_frac=BASELINE_EDGE_FRACTION,
    manual_pW=MANUAL_BASELINE_PW_PER_HZ
)
y_above_raw_pW = np.clip(y_win_raw_pW - baseline_raw_pW, 0.0, None)  # only above the line
area_raw_W = integrate_trapz(x_win_raw_MHz * 1e6, y_above_raw_pW * 1e-12)
area_raw_pW = area_raw_W * 1e12

# --------------- averaged integration (same baseline logic) ---------------
area_avg_W = None
area_avg_pW = None
x_win_avg_MHz = None
y_above_avg_pW = None
if x_avg.size and y_avg.size:
    mask_avg = (x_avg >= f0) & (x_avg <= f1)
    if np.any(mask_avg):
        x_win_avg_MHz = x_avg[mask_avg]
        y_win_avg_pW = y_avg[mask_avg]
        baseline_avg_pW = estimate_baseline(
            y_win_avg_pW, mode=BASELINE_MODE.lower(),
            edge_frac=BASELINE_EDGE_FRACTION,
            manual_pW=MANUAL_BASELINE_PW_PER_HZ
        )
        y_above_avg_pW = np.clip(y_win_avg_pW - baseline_avg_pW, 0.0, None)
        area_avg_W = integrate_trapz(x_win_avg_MHz * 1e6, y_above_avg_pW * 1e-12)
        area_avg_pW = area_avg_W * 1e12

print("\n=== Peak Area (above baseline) ===")
print(f"Window: [{f0:.6f}, {f1:.6f}] MHz  | Baseline mode: {BASELINE_MODE}")
if BASELINE_MODE.lower() == 'manual':
    print(f"Manual baseline: {MANUAL_BASELINE_PW_PER_HZ:.3g} pW/Hz")
print(f"RAW      area = {area_raw_W:.6e} W   = {area_raw_pW:.6e} pW")
if area_avg_W is not None:
    print(f"AVERAGED area = {area_avg_W:.6e} W   = {area_avg_pW:.6e} pW")

# --------------- plot ---------------
fig, ax = plt.subplots(figsize=(10, 6))

# averaged curve
if x_avg.size and y_avg.size:
    ax.plot(x_avg, y_avg, '-o', ms=3.5, lw=1.8, label=f"{AVG_BLOCK}-point average")

# show manual baseline as dotted line
if BASELINE_MODE.lower() == 'manual':
    ax.axhline(MANUAL_BASELINE_PW_PER_HZ, color='k', ls=':', lw=1.2, alpha=0.9,
               label=f"manual baseline = {MANUAL_BASELINE_PW_PER_HZ:g} pW/Hz")

# --- RED RAW FILL (still disabled) ---
# To visualize the RAW area above baseline in RED, uncomment the block below.
# x_fill_raw = x_win_raw_MHz.copy()
# y_fill_raw_pW = y_above_raw_pW.copy()  # already baseline-subtracted and clipped at 0
# order = np.argsort(x_fill_raw)
# x_fill_raw = x_fill_raw[order]
# y_fill_raw_pW = y_fill_raw_pW[order]
# ax.fill_between(x_fill_raw, MANUAL_BASELINE_PW_PER_HZ if BASELINE_MODE.lower()=='manual' else 0.0,
#                 (y_fill_raw_pW + (MANUAL_BASELINE_PW_PER_HZ if BASELINE_MODE.lower()=='manual' else 0.0)),
#                 color='crimson', alpha=0.35, label='Integrated area (RAW, red)')

# --- GREEN AVERAGED FILL (shown) ---
if FILL_AVG_GREEN and (x_win_avg_MHz is not None) and (y_above_avg_pW is not None):
    order = np.argsort(x_win_avg_MHz)
    x_fill_avg = x_win_avg_MHz[order]
    y_fill_avg_pW = y_above_avg_pW[order]
    y_base = MANUAL_BASELINE_PW_PER_HZ if BASELINE_MODE.lower() == 'manual' else 0.0
    ax.fill_between(x_fill_avg, y_base, y_base + y_fill_avg_pW,
                    color='green', alpha=0.30, label='Integrated area (AVG, green)')

ax.set_xlabel("Frequency (MHz)", fontsize=18, fontname="Calibri")
ax.set_ylabel("PSD (pW/Hz)", fontsize=18, fontname="Calibri")
ax.set_title("Noise-subtracted PSD: area above a horizontal baseline", fontsize=18, fontname="Calibri")

ax.tick_params(axis='both', which='major', labelsize=16)
for lab in (ax.get_xticklabels() + ax.get_yticklabels()):
    lab.set_fontname("Calibri")

annot = f"RAW area = {area_raw_pW:.3e} pW"
if area_avg_pW is not None:
    annot += f"\nAVG area = {area_avg_pW:.3e} pW"
annot += f"\n[{PEAK_START_MHZ:.6f}, {PEAK_END_MHZ:.6f}] MHz | mode: {BASELINE_MODE}"
if BASELINE_MODE.lower() == 'manual':
    annot += f"\nmanual baseline: {MANUAL_BASELINE_PW_PER_HZ:g} pW/Hz"
ax.text(0.02, 0.98, annot, transform=ax.transAxes, ha='left', va='top')

ax.legend(frameon=False)
out_png = os.path.join(out_dir, f"psd_area_above_baseline_run{run_data}_minus_{run_background}_avg{AVG_BLOCK}.png")
plt.savefig(out_png, dpi=220, bbox_inches='tight')
print("Saved figure to:", out_png)

try:
    plt.show()
except Exception as e:
    print("plt.show() failed:", e)



