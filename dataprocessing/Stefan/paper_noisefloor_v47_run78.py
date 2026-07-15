import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
from database import *

# v47 baseband spectra, full compensated traces:
# fig 1 = run 78 (true centre 0); figs 2+3 = runs 78+82+84 concatenated (log-x and linear-x)

F_MIX, RBW, WIN, GLITCH_HW = 2.5e6, 1.676, 1001, 1   # WIN: broad rolling avg -> smooth filter shape
FMIN = 1.0
GAIN = 200 * 5.65   # amplification chain
Z    = 7.5e3        # Ohm, V -> I conversion
DB, REF_RUN = r'D:\databases CD12_B5_F4\CD12_B5_F4v47_08_12_25.db', 72
CONCAT_RUNS = [78, 82, 84]   # true centres 0, 100 kHz, 150 kHz

def get_spec(run_id):
    f, s = extract_1d(run_id, data_1d_name="V_fft_avg_avg",
                      setpoint_name="freq_param", plot=False, return_exp_name=False)
    return np.asarray(f, float), np.asarray(s, float)

def rolling_mean(y, win):
    pad = win // 2
    return np.convolve(np.pad(y, pad, mode="edge"), np.ones(win) / win, mode="valid")

qc.config["core"]["db_location"] = DB
fref, sref = get_spec(REF_RUN)
clean = sref.copy()
i0 = np.argmax(clean)
clean[max(0, i0 - GLITCH_HW):i0 + GLITCH_HW + 1] = np.nan
bad = np.isnan(clean)
clean[bad] = np.interp(fref[bad], fref[~bad], clean[~bad])
ref = rolling_mean(clean, WIN)
ref_norm = ref / ref.max()

def spec_comp(run_id):
    f, s = get_spec(run_id)
    R = ref_norm if s.size == ref_norm.size else np.interp(
        np.linspace(0, 1, s.size), np.linspace(0, 1, ref_norm.size), ref_norm)
    ft = f - F_MIX
    m = ft >= FMIN
    return ft[m], (s[m] / R[m]) / np.sqrt(RBW) / (GAIN * Z) * 1e15   # fA/sqrt(Hz)

plt.rcParams.update({
    "font.size": 8, "axes.linewidth": 0.6, "legend.frameon": False,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})

def make_fig(x, y, label, fname, xscale="log", xlabel="frequency (Hz)"):
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax.plot(x, y, color="#1baf7a", lw=0.4, alpha=0.9, label=label)
    ax.set_xscale(xscale); ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\sqrt{S_I}$ (fA/$\sqrt{\mathrm{Hz}}$)")
    ax.grid(True, which="major", color="#e1e0d9", lw=0.4)
    ax.set_axisbelow(True)
    ax.legend(fontsize=7, loc="upper right", handletextpad=0.4)
    plt.tight_layout(pad=0.3)
    plt.savefig(fname + ".png", dpi=300, bbox_inches="tight")
    plt.savefig(fname + ".pdf", bbox_inches="tight")

import os
FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)
BASE = os.path.join(FIGDIR, os.path.basename(__file__).replace(".py", ""))
ftc, asdc = spec_comp(CONCAT_RUNS[0])
make_fig(ftc, asdc, "v47 run 78 (0 MHz)", BASE)

for run_id in CONCAT_RUNS[1:]:            # concatenate: each run contributes above the previous edge
    ft, asd = spec_comp(run_id)
    m = ft > ftc.max()
    ftc, asdc = np.concatenate([ftc, ft[m]]), np.concatenate([asdc, asd[m]])
lab = "v47 runs " + "+".join(str(r) for r in CONCAT_RUNS)
make_fig(ftc, asdc, lab, BASE + "_concat")
make_fig(ftc / 1e3, asdc, lab, BASE + "_concat_linx",
         xscale="linear", xlabel="frequency (kHz)")
plt.show()
