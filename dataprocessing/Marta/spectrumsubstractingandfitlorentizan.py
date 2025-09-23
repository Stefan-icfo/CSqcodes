#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit

# ===================== USER CONFIG =====================
run_data = 1480
run_background = 1465
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v8.db"
out_dir  = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Figures"
os.makedirs(out_dir, exist_ok=True)

# Subtraction mode:
DO_VOLTAGE_DOMAIN = True   # True: V_sig = sqrt(P_data) - sqrt(P_bg);  P_sig = V_sig^2
USE_INDEX_AXIS    = False  # True: plot x = point index; False: plot x = data-run frequency (truncated)

# Optional: restrict the fit to a window (on the x axis you’re using); set None to use all.
FIT_WINDOW = None   # e.g. (141.85, 141.97) if using frequency; or (200, 400) if using index
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

# Lorentzian model
def lorentzian(x, x0, gamma, A, y0):
    # y = y0 + A*gamma^2 / ((x-x0)^2 + gamma^2)
    return y0 + (A * gamma**2) / ((x - x0)**2 + gamma**2)

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
    y_unit = "W/Hz"
else:
    y_sig = y_d[:n] - y_b[:n]
    y_unit = "W/Hz"

# Convert to pW/Hz for your plot if you like
y_sig_plot = y_sig / 1e-12   # pW/Hz
y_label = f"PSD (pW/Hz)"

print(f"[subtract] using first {n} points (index-by-index).")
print(f"  result rng: [{np.nanmin(y_sig):.3e}, {np.nanmax(y_sig):.3e}] {y_unit}")

# --------------- Lorentzian FIT on subtracted data ---------------
# Choose fit data (optionally windowed)
mask_fit = np.isfinite(x_plot) & np.isfinite(y_sig_plot)
if FIT_WINDOW is not None:
    lo, hi = FIT_WINDOW
    mask_fit &= (x_plot >= lo) & (x_plot <= hi)

xf = x_plot[mask_fit]
yf = y_sig_plot[mask_fit]

fit_ok = False
popt = perr = None

if xf.size >= 8 and np.nanmax(yf) > np.nanmin(yf):
    # Initial guesses
    x0_guess = xf[np.nanargmax(yf)]
    span = max(1e-9, (xf.max() - xf.min()))
    gamma_guess = span / 50.0           # broad but reasonable default
    A_guess = float(np.nanmax(yf) - np.nanmedian(yf))
    y0_guess = float(np.nanmedian(yf))
    p0 = [x0_guess, gamma_guess, A_guess, y0_guess]

    # Bounds: gamma >= 0
    bounds = (
        [xf.min(),     0.0, -np.inf, -np.inf],
        [xf.max(),  np.inf,  np.inf,  np.inf],
    )

    try:
        popt, pcov = curve_fit(lorentzian, xf, yf, p0=p0, bounds=bounds, maxfev=20000)
        perr = np.sqrt(np.diag(pcov))
        fit_ok = True
        print("\n✅ Lorentzian fit (on subtracted data, plotted units):")
        print(f"  x0     = {popt[0]:.6f} ± {perr[0]:.2e}")
        print(f"  gamma  = {popt[1]:.6f} ± {perr[1]:.2e}")
        print(f"  A      = {popt[2]:.3e} ± {perr[2]:.1e}")
        print(f"  y0     = {popt[3]:.3e} ± {perr[3]:.1e}")
    except Exception as e:
        print(f"⚠️  Fit failed: {e}")
else:
    print("⚠️  Not enough valid points for fitting; skipping fit.")

# --------------- plot ---------------
plt.figure(figsize=(10, 6))
plt.plot(x_plot, y_sig_plot, '-', lw=1.8,
         label=f"Corrected: run {run_data} − {run_background}")

# overlay fit if available
if fit_ok:
    xx = np.linspace(xf.min(), xf.max(), 1500)
    yy = lorentzian(xx, *popt)
    plt.plot(xx, yy, '-', lw=2.2, color='crimson', label='Lorentzian fit')

xlabel = "Point index" if USE_INDEX_AXIS or fkey_d is None else "Frequency (MHz)"
plt.xlabel(xlabel, fontsize=18, fontname="Calibri")
plt.ylabel(y_label, fontsize=18, fontname="Calibri")
title = "Point-by-point subtraction (ignoring frequency mismatch)"
if FIT_WINDOW is not None:
    title += f" | fit in [{FIT_WINDOW[0]}, {FIT_WINDOW[1]}]"
plt.title(title, fontsize=18, fontname="Calibri")

plt.tight_layout()

plt.tick_params(axis='both', which='major', labelsize=16)
for lab in (plt.gca().get_xticklabels() + plt.gca().get_yticklabels()):
    lab.set_fontname("Calibri")

out_png = os.path.join(out_dir, f"psd_point_by_point_run{run_data}_minus_{run_background}_withfit.png")
plt.savefig(out_png, dpi=200)
print("Saved figure to:", out_png)

try:
    plt.show()
except Exception as e:
    print("plt.show() failed:", e)
