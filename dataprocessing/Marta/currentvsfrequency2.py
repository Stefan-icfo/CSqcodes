#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot CURRENT vs FREQUENCY for run 873 from QCoDeS database.
- Use the parameters defined in mech_simple_fun_db: frequency sweep param and 'I_rf'.
- Plot ALL raw data points (current vs frequency).
- Additionally, compute block averages (every 10 points) and plot them connected by a red line.

Author: @marta.cagetti
"""

import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from pathlib import Path

# =============================
# User settings
# =============================
DB_PATH = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v7.db"
RUN_ID  = 909

CURR_KEY = 'I_rf'       # current registered in mech_simple_fun_db
BLOCK_SIZE = 60       # average every 10 points

# =============================
# Helpers
# =============================

def _flatten(a):
    try:
        return np.asarray(a).astype(float).ravel()
    except Exception:
        return None

def get_fx_I_from_dataset(ds, curr_key=CURR_KEY):
    """
    Extract frequency setpoint (x) and I_rf (y) from dataset using parameter_data.
    """
    pd = ds.get_parameter_data()
    if curr_key in pd:
        blk = pd[curr_key]
        # find the setpoint key
        sp_keys = [k for k in blk.keys() if k != curr_key]
        if sp_keys:
            x = _flatten(blk[sp_keys[0]])
            y = _flatten(blk[curr_key])
            return x, y
    raise KeyError(f'Could not find {curr_key} in parameter_data.')

# =============================
# Load and process
# =============================
qc.config["core"]["db_location"] = DB_PATH

ds = qc.load_by_id(RUN_ID)
print(f"[QCoDeS] Loaded run {RUN_ID}.")

x, y = get_fx_I_from_dataset(ds)
mask = np.isfinite(x) & np.isfinite(y)
x, y = x[mask], y[mask]

# =============================
# Compute block averages of 10 points
# =============================
N = len(x)
n_blocks = N // BLOCK_SIZE
xb, yb = [], []
for i in range(n_blocks):
    s, e = i*BLOCK_SIZE, (i+1)*BLOCK_SIZE
    xb.append(np.mean(x[s:e]))
    yb.append(np.mean(y[s:e]))
xb = np.array(xb)
yb = np.array(yb)

# =============================
# Plot RAW DATA and AVERAGED BLOCKS
# =============================
x_plot = x.copy()
xb_plot = xb.copy()
x_label = 'Frequency [Hz]'
if np.nanmedian(np.abs(x_plot)) > 1e6:
    x_plot = x_plot / 1e6
    xb_plot = xb_plot / 1e6
    x_label = 'Frequency [MHz]'

# Helper: average in consecutive chunks of size N
def chunk_average(xa, ya, n=60):
    m = len(xa) // n
    if m == 0:
        return np.array([]), np.array([])
    xa = xa[:m*n]
    ya = ya[:m*n]
    xb = xa.reshape(m, n).mean(axis=1)
    yb = ya.reshape(m, n).mean(axis=1)
    return xb, yb

# Sort by frequency for visual clarity
idx = np.argsort(x)
x_sorted = x[idx]
y_sorted = y[idx]

# Compute 10-point chunk averages on the sorted data
x_avg, y_avg = chunk_average(x_sorted, y_sorted, n=60)

# Convert to MHz if clearly RF
x_plot = x_sorted.copy()
x_label = 'Frequency [Hz]'
if np.nanmedian(np.abs(x_plot)) > 1e6:
    x_plot = x_plot / 1e6
    x_avg_plot = x_avg / 1e6
    x_label = 'Frequency [MHz]'
else:
    x_avg_plot = x_avg

plt.figure(figsize=(8.5,5.2))
# raw data (blue)
plt.plot(x_plot, y_sorted, 'o', ms=2, alpha=0.6, label='raw')
# 10-point averaged (red line with markers)
if x_avg_plot.size:
    plt.plot(x_avg_plot, y_avg, '-o', linewidth=1.8, markersize=3.5, color='red', label='10-point avg')

plt.xlabel(x_label)
plt.ylabel('Current [A]')
plt.title(f'Run {RUN_ID}: I_rf vs frequency (raw + pt avg)')
plt.grid(True, ls=':', alpha=0.6)
plt.legend()
plt.tight_layout()

out = Path.cwd() / f"run_{RUN_ID}_Irf_vs_f_with_pt_avg.png"
plt.savefig(out, dpi=200)
plt.show()
print("Saved:", out)






