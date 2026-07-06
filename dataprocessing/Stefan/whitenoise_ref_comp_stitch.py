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
GLITCH_HW = 1        # half-width of bins removed around the outlier spike
RBW       = 1.676    # Hz, resolution bandwidth -> calibrates amplitude to V/sqrt(Hz)
STEP      = 50e3     # Hz, run-to-run centre spacing = width of each stitched tile
RUNS      = [152, 155, 158,161]   # runs to plot / stitch (edit to taste)

ASD_LABEL = r"ASD [V/$\sqrt{\mathrm{Hz}}$]"

USE_MIDPOINT_AS_CENTER = True     # True if freq_param is absolute; False to read centre from metadata
CENTER_FIELD = "center_freq"

def get_spec(run_id):
    f, s = extract_1d(run_id, data_1d_name="V_fft_avg_avg",
                      setpoint_name="freq_param", plot=False, return_exp_name=False)
    return np.asarray(f, float), np.asarray(s, float)

def get_center(run_id, f):
    if USE_MIDPOINT_AS_CENTER:
        return 0.5 * (f.min() + f.max())
    return get_metadata(run_id - 1, print_it=False, return_data=True)[CENTER_FIELD]

def rolling_mean(y, win):
    if win % 2 == 0:
        win += 1
    pad = win // 2
    return np.convolve(np.pad(y, pad, mode="edge"), np.ones(win) / win, mode="valid")

# ---- build the reference shape from run 188 ----
fref, sref = get_spec(REF_RUN)
clean = sref.copy()
i0 = np.argmax(clean)
lo, hi = max(0, i0 - GLITCH_HW), min(clean.size, i0 + GLITCH_HW + 1)
clean[lo:hi] = np.nan
bad = np.isnan(clean)
clean[bad] = np.interp(fref[bad], fref[~bad], clean[~bad])
ref_smooth = rolling_mean(clean, WIN)
ref_norm   = ref_smooth / ref_smooth.max()

# ---- compensate one run, return calibrated ASD ----
def compensate(run_id, show=True):
    f, s = get_spec(run_id)
    if s.size == ref_norm.size:
        R = ref_norm
    else:
        R = np.interp(np.linspace(0, 1, s.size),
                      np.linspace(0, 1, ref_norm.size), ref_norm)
    asd     = (s / R) / np.sqrt(RBW)
    asd_raw =  s      / np.sqrt(RBW)

    if show:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(f, asd_raw, lw=0.6, alpha=0.6, label=f"run {run_id} raw")
        ax.plot(f, asd,     lw=0.8,            label=f"run {run_id} compensated")
        ax.set_yscale("log")
        ax.set_xlabel("freq_param [Hz]"); ax.set_ylabel(ASD_LABEL)
        ax.set_title(f"run {run_id} flattened by run-{REF_RUN} shape"); ax.legend()
        plt.tight_layout(); plt.show()

        off = f - get_center(run_id, f)
        m = off <= 0
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(off[m], asd_raw[m], lw=0.6, alpha=0.6, label="raw")
        ax.plot(off[m], asd[m],     lw=0.8,            label="compensated")
        ax.set_yscale("log")
        ax.set_xlabel("frequency offset from centre [Hz]"); ax.set_ylabel(ASD_LABEL)
        ax.set_title(f"run {run_id} — lower half (centre = 0)"); ax.legend()
        plt.tight_layout(); plt.show()

    return f, asd

# ---- reference check ----
plt.figure(figsize=(7, 3))
plt.plot(fref, sref, lw=0.5, alpha=0.5, label=f"run {REF_RUN} raw")
plt.plot(fref, ref_smooth, "r", lw=1.3, label="smoothed reference (deglitched)")
plt.legend(); plt.title("white-noise reference shape")
plt.xlabel("freq_param [Hz]"); plt.tight_layout(); plt.show()

# ---- individual runs ----
for run_id in RUNS:
    compensate(run_id, show=True)

# ---- last plot: stitched composite, first run's lower half starting at 0 ----
f0, _ = get_spec(RUNS[0])
ref0 = get_center(RUNS[0], f0) - STEP        # lower EDGE of first tile -> x = 0
fig, ax = plt.subplots(figsize=(8, 4))
for run_id in RUNS:
    f, asd = compensate(run_id, show=False)
    c = get_center(run_id, f)
    off = f - c
    seg = (off <= 0) & (off >= -STEP)         # lower half, one tile wide
    ax.plot(f[seg] - ref0, asd[seg], lw=0.7,
            label=f"run {run_id}  ({c/1e6:.3f} MHz)")
ax.set_yscale("log")
ax.set_xlabel(f"stitched frequency [Hz]  (run {RUNS[0]} lower edge = 0)")
ax.set_ylabel(ASD_LABEL)
ax.set_title("stitched spectrum (lower-half tiles)")
ax.legend(); plt.tight_layout(); plt.show()