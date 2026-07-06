import numpy as np
import matplotlib.pyplot as plt
import scipy as scp
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

qc.config["core"]["db_location"] = r'D:\databases CD12_B5_F4\CD12_B5_F4v49_15_12_25.db'

REF_RUN   = 188      # white-noise reference: flat input -> pure system response
WIN       = 51       # rolling-average window (odd); larger = smoother
GLITCH_HW = 1        # half-width of bins removed around the outlier spike (1 -> 3 bins)
RBW       = 1.676    # Hz, resolution bandwidth -> calibrates amplitude to V/sqrt(Hz)

ASD_LABEL = r"ASD [V/$\sqrt{\mathrm{Hz}}$]"

def get_spec(run_id):
    f, s = extract_1d(run_id, data_1d_name="V_fft_avg_avg",
                      setpoint_name="freq_param", plot=False, return_exp_name=False)
    return np.asarray(f, float), np.asarray(s, float)

def rolling_mean(y, win):
    if win % 2 == 0:
        win += 1
    pad = win // 2
    return np.convolve(np.pad(y, pad, mode="edge"), np.ones(win) / win, mode="valid")

# ---- build the reference shape from run 188 ----
fref, sref = get_spec(REF_RUN)
clean = sref.copy()
i0 = np.argmax(clean)                                   # the one high outlier
lo, hi = max(0, i0 - GLITCH_HW), min(clean.size, i0 + GLITCH_HW + 1)
clean[lo:hi] = np.nan                                   # remove it
bad = np.isnan(clean)
clean[bad] = np.interp(fref[bad], fref[~bad], clean[~bad])   # bridge the gap

ref_smooth = rolling_mean(clean, WIN)                   # shape-preserving smooth
ref_norm   = ref_smooth / ref_smooth.max()              # divide targets by this

# ---- reusable compensator: target spectrum * (1 / reference shape), in V/sqrt(Hz) ----
def compensate(run_id, show=True):
    f, s = get_spec(run_id)
    if s.size == ref_norm.size:
        R = ref_norm
    else:
        # align by position across the band (both centred demod FFTs).
        # If the runs share an ABSOLUTE-Hz grid instead, use:  R = np.interp(f, fref, ref_norm)
        R = np.interp(np.linspace(0, 1, s.size),
                      np.linspace(0, 1, ref_norm.size), ref_norm)
    comp     = s / R
    asd      = comp / np.sqrt(RBW)        # V/sqrt(Hz)
    asd_raw  = s    / np.sqrt(RBW)

    if show:
        # --- plot 1: full band ---
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(f, asd_raw, lw=0.6, alpha=0.6, label=f"run {run_id} raw")
        ax.plot(f, asd,     lw=0.8,            label=f"run {run_id} compensated")
        ax.set_yscale("log")
        ax.set_xlabel("freq_param [Hz]"); ax.set_ylabel(ASD_LABEL)
        ax.set_title(f"run {run_id} flattened by run-{REF_RUN} white-noise shape")
        ax.legend(); plt.tight_layout(); plt.show()

        # --- plot 2: upper half only, band centre re-referenced to 0 ---
        fc  = 0.5 * (f.min() + f.max())
        off = f - fc
        m   = off >= 0
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(off[m], asd_raw[m], lw=0.6, alpha=0.6, label="raw")
        ax.plot(off[m], asd[m],     lw=0.8,            label="compensated")
        ax.set_yscale("log")
        ax.set_xlabel("frequency offset from centre [Hz]"); ax.set_ylabel(ASD_LABEL)
        ax.set_title(f"run {run_id} — upper half (centre = 0)")
        ax.legend(); plt.tight_layout(); plt.show()

    return f, comp

# ---- check the reference itself ----
plt.figure(figsize=(7, 3))
plt.plot(fref, sref, lw=0.5, alpha=0.5, label=f"run {REF_RUN} raw")
plt.plot(fref, ref_smooth, "r", lw=1.3, label="smoothed reference (deglitched)")
plt.legend(); plt.title("white-noise reference shape")
plt.xlabel("freq_param [Hz]"); plt.tight_layout(); plt.show()

# ---- start with run 152 ----
compensate(152)