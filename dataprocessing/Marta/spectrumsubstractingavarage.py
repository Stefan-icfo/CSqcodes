#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ===================== USER CONFIG =====================
run_data = 3971
run_background = 3977
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v11.db"
out_dir  = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Figures"
os.makedirs(out_dir, exist_ok=True)

# Subtraction mode:
DO_VOLTAGE_DOMAIN = True   # True: V_sig = sqrt(P_data) - sqrt(P_bg);  P_sig = V_sig^2
USE_INDEX_AXIS    = False  # True: plot x = point index; False: plot x = data-run frequency (truncated)

# Averaging:
AVG_BLOCK = 5    # average each 8 consecutive points (non-overlapping)
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
    """
    Load a 1D trace in ACQUISITION ORDER (no sorting).
    Returns (x_MHz, y, dep_key, freq_key).
    """
    ds = qc.load_by_id(run_id)
    pd = ds.get_parameter_data()

    # choose dependent
    dep = next((k for k in CANDIDATE_SIGNAL_KEYS if k in pd), next(iter(pd.keys())))
    blk = pd[dep]

    # choose frequency setpoint inside the same block (only for plotting x)
    freq_key = None
    for k in blk.keys():
        if k == dep:
            continue
        if any(t in k.lower() for t in FREQ_KEY_CANDIDATES):
            freq_key = k
            break
    if freq_key is None:
        # just pick any other setpoint if frequency key name is unusual
        freq_key = next((k for k in blk.keys() if k != dep), None)

    # flatten in given order (no sort)
    x = np.asarray(blk[freq_key]).ravel() / 1e6 if freq_key else np.arange(len(np.asarray(blk[dep]).ravel()), dtype=float)
    y = np.asarray(blk[dep]).ravel()

    # clean NaNs (keep index alignment)
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
    """
    Non-overlapping block average of size n.
    Tail points that don't fill a full block are dropped.
    """
    if n <= 1:
        return x.copy(), y.copy()
    m = (len(y) // n) * n
    if m == 0:
        return np.array([]), np.array([])
    xb = x[:m].reshape(-1, n).mean(axis=1)
    yb = y[:m].reshape(-1, n).mean(axis=1)
    return xb, yb

# --------------- load both runs ---------------
x_d, y_d, dep_d, fkey_d = load_1d_in_order(run_data)
x_b, y_b, dep_b, fkey_b = load_1d_in_order(run_background)

# --------------- point-by-point subtraction (index-based) ---------------
n = min(len(y_d), len(y_b))
if n == 0:
    raise RuntimeError("No points to subtract (one of the traces is empty).")

x_plot = (np.arange(n, dtype=float) if USE_INDEX_AXIS else x_d[:n])

if DO_VOLTAGE_DOMAIN:
    V_d = np.sqrt(np.clip(y_d[:n], 0.0, None))
    V_b = np.sqrt(np.clip(y_b[:n], 0.0, None))
    V_sig = V_d - V_b
    y_sig = V_sig**2
    y_label = "PSD (pW/Hz)"
else:
    y_sig = y_d[:n] - y_b[:n]
    y_label = "PSD (W/Hz)  [power subtraction]"

print(f"[subtract] using first {n} points (index-by-index).")
print(f"  data range: [{np.nanmin(y_d[:n]):.3e}, {np.nanmax(y_d[:n]):.3e}]")
print(f"  back range: [{np.nanmin(y_b[:n]):.3e}, {np.nanmax(y_b[:n]):.3e}]")
print(f"  result rng: [{np.nanmin(y_sig):.3e}, {np.nanmax(y_sig):.3e}]")

# Convert to pW/Hz for plotting (your convention)
y_sig_plot = y_sig / 1e-12

# --------------- 8-point averaging ---------------
x_avg, y_avg = block_average(x_plot, y_sig_plot, AVG_BLOCK)

# --------------- plot ---------------
plt.figure(figsize=(10, 6))

# raw corrected spectrum# 
# plt.plot(x_plot, y_sig_plot, '-', lw=1.2, alpha=0.6,
        #  label=f"Corrected: run {run_data} − {run_background}")

# averaged curve
if x_avg.size and y_avg.size:
    plt.plot(x_avg, y_avg, '-o', ms=3.5, lw=1.8,
             label=f"{AVG_BLOCK}-point average")

xlabel = "Point index" if USE_INDEX_AXIS or fkey_d is None else "Frequency (MHz)"
plt.xlabel(xlabel, fontsize=18, fontname="Calibri")
plt.ylabel(y_label, fontsize=18, fontname="Calibri")
plt.title("Point-by-point subtraction with 8-point averaging", fontsize=18, fontname="Calibri")



plt.tick_params(axis='both', which='major', labelsize=16)
for lab in (plt.gca().get_xticklabels() + plt.gca().get_yticklabels()):
    lab.set_fontname("Calibri")

out_png = os.path.join(out_dir, f"psd_point_by_point_run{run_data}_minus_{run_background}_avg{AVG_BLOCK}.png")
plt.savefig(out_png, dpi=200)
print("Saved figure to:", out_png)

try:
    plt.show()
except Exception as e:
    print("plt.show() failed:", e)


