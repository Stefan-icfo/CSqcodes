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

qc.config["core"]["db_location"] = r'D:\databases CD12_B5_F4\CD12_B5_F4v5.db'

# v5 (Aug 2025), older notation: freq_param is absolute (no 2.5 MHz mix-down).
# Weak-drive frequency sweep across the ~138 MHz mode ("thermal spectrum+mix_and_weakly_driven_at_X"):
# the drive tone steps 5 kHz per run and the 109.86 kHz window follows it, so each run holds
# the thermal spectrum + the driven peak at window centre.
# No usable white-noise/background run with matching filter settings exists (run 156's
# 220 kHz background has a different relative filter shape), so the window shape is
# self-referenced: median across all runs per bin index (features move run-to-run,
# the common window response stays).
WIN       = 51       # rolling-average window (odd); larger = smoother
GLITCH_HW = 1        # half-width of bins removed around the outlier spike
RBW       = 1.676    # Hz, resolution bandwidth (= tile bin spacing) -> calibrates amplitude to V/sqrt(Hz)
RUNS_A    = [274, 276, 278, 280, 282, 284, 286, 288, 292, 294, 296,
             298, 304, 306, 308, 310, 312, 314, 316]   # sweep A: drive 137.950..138.050 MHz, 5 kHz steps
                                                       # (300+302 excluded: glitch, floor drops to 0.55x/0.90x)
RUNS_B    = [318, 320, 322, 324, 326, 328, 330, 332, 334, 336,
             338, 340, 342, 344, 346, 348, 350, 352, 354]   # sweep B repeat: 137.950..138.040 MHz; floor 1.21x sweep A (drift)
RUNS      = RUNS_A + RUNS_B
SHOW_INDIVIDUAL = False   # per-run before/after plots; False avoids 80 popups for the full set

ASD_LABEL = r"ASD [V/$\sqrt{\mathrm{Hz}}$]"

def get_spec(run_id):
    f, s = extract_1d(run_id, data_1d_name="V_fft_avg_avg",
                      setpoint_name="freq_param", plot=False, return_exp_name=False)
    return np.asarray(f, float), np.asarray(s, float)

def get_center(f):
    return 0.5 * (f.min() + f.max())

def rolling_mean(y, win):
    if win % 2 == 0:
        win += 1
    pad = win // 2
    return np.convolve(np.pad(y, pad, mode="edge"), np.ones(win) / win, mode="valid")

# ---- build the self-referenced window shape: median across tiles per bin index ----
SPECS = {run_id: get_spec(run_id) for run_id in RUNS}
stack = np.vstack([s for _, s in SPECS.values()])
med   = np.median(stack, axis=0)
i0 = np.argmax(med)                       # drive tone sits at window centre in every tile
lo, hi = max(0, i0 - GLITCH_HW), min(med.size, i0 + GLITCH_HW + 1)
med[lo:hi] = np.nan
bad = np.isnan(med)
med[bad] = np.interp(np.flatnonzero(bad), np.flatnonzero(~bad), med[~bad])
ref_smooth = rolling_mean(med, WIN)
ref_norm   = ref_smooth / ref_smooth.max()

# ---- compensate one run, return calibrated ASD ----
def compensate(run_id, show=True):
    f, s = SPECS[run_id]
    R = ref_norm if s.size == ref_norm.size else np.interp(
        np.linspace(0, 1, s.size), np.linspace(0, 1, ref_norm.size), ref_norm)
    asd     = (s / R) / np.sqrt(RBW)
    asd_raw =  s      / np.sqrt(RBW)

    if show:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(f, asd_raw, lw=0.6, alpha=0.6, label=f"run {run_id} raw")
        ax.plot(f, asd,     lw=0.8,            label=f"run {run_id} compensated")
        ax.set_yscale("log")
        ax.set_xlabel("freq_param [Hz]"); ax.set_ylabel(ASD_LABEL)
        ax.set_title(f"run {run_id} flattened by self-referenced shape"); ax.legend()
        plt.tight_layout(); plt.show()

    return f, asd

# ---- reference check ----
plt.figure(figsize=(7, 3))
plt.plot(np.median(stack, axis=0), lw=0.5, alpha=0.5, label="median across tiles (raw)")
plt.plot(ref_smooth, "r", lw=1.3, label="smoothed reference (deglitched)")
plt.legend(); plt.title("self-referenced window shape (median over tiles)")
plt.xlabel("bin index"); plt.tight_layout(); plt.show()

# ---- individual runs ----
for run_id in RUNS:
    compensate(run_id, show=SHOW_INDIVIDUAL)

# ---- composite: all windows at absolute frequency, coloured by drive frequency ----
from matplotlib import cm
from matplotlib.colors import Normalize
F0 = 138e6
centers = {r: get_center(SPECS[r][0]) for r in RUNS}
off_k = {r: (c - F0) / 1e3 for r, c in centers.items()}
cnorm = Normalize(min(off_k.values()), max(off_k.values()))
fig, ax = plt.subplots(figsize=(9, 4.5))
for run_id in RUNS:
    f, asd = compensate(run_id, show=False)
    ax.plot(f / 1e6, asd, lw=0.6, alpha=0.8,
            color=cm.viridis(cnorm(off_k[run_id])),
            ls="-" if run_id in RUNS_A else "--")
ax.set_yscale("log")
ax.set_xlabel("frequency [MHz]"); ax.set_ylabel(ASD_LABEL)
ax.set_title(f"weak-drive sweep across the ~138 MHz mode (v5, {len(RUNS)} runs)")
fig.colorbar(cm.ScalarMappable(norm=cnorm, cmap="viridis"), ax=ax,
             label="drive offset from 138 MHz [kHz]")
ax.plot([], [], "k-", label="sweep A"); ax.plot([], [], "k--", label="sweep B")
ax.legend(fontsize=8); plt.tight_layout(); plt.show()
