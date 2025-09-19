#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot CURRENT vs FREQUENCY for run 873 from QCoDeS database.
- Use the parameters defined in mech_simple_fun_db: frequency sweep param and 'I_rf'.
- Plot ALL raw data points (current vs frequency).

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
RUN_ID  = 873

CURR_KEY = 'I_rf'       # current registered in mech_simple_fun_db

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
# Plot RAW DATA
# =============================
x_plot = x.copy()
x_label = 'Frequency [Hz]'
if np.nanmedian(np.abs(x_plot)) > 1e6:
    x_plot = x_plot / 1e6
    x_label = 'Frequency [MHz]'

plt.figure(figsize=(8,5))
plt.plot(x_plot, y, 'o', ms=2, alpha=0.7)
plt.xlabel(x_label)
plt.ylabel('Current [A]')
plt.title(f'Run {RUN_ID}: I_rf vs frequency')
plt.grid(True, ls=':', alpha=0.6)
plt.tight_layout()

out = Path.cwd() / f"run_{RUN_ID}_Irf_vs_f.png"
plt.savefig(out, dpi=200)
plt.show()
print("Saved:", out)




