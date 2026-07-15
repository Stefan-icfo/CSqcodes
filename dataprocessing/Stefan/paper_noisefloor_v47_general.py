import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
from database import *

# v47 baseband spectrum, fully concatenated compensated trace: splice consecutive
# tiles (each contributes above the previous edge) until the first coverage hole.

F_MIX, RBW, WIN, GLITCH_HW = 2.5e6, 1.676, 1001, 1   # WIN: broad rolling avg -> smooth filter shape
FMIN = 1.0
GAIN = 200 * 5.65   # amplification chain
Z    = 7.5e3        # Ohm, V -> I conversion
PROBE_F = 100e3     # Hz, probe tone: remove the spike at this frequency (peak +- GLITCH_HW bins)
DB, REF_RUN = r'D:\databases CD12_B5_F4\CD12_B5_F4v47_08_12_25.db', 72
CANDIDATES = [78, 80, 82, 84, 86, 88, 90, 93, 95, 97, 99, 101, 103, 105, 107, 109,
              111, 114, 116, 118,                      # linear ladder: 0..950 kHz, 50 kHz steps
              9, 11, 13, 15, 17, 19, 21, 23, 25, 27]   # log ladder: 1 MHz.. (holes start here)

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

specs, top = [], None
for run_id in CANDIDATES:
    ft, asd = spec_comp(run_id)
    if top is not None and ft.min() > top:
        print(f"hole before run {run_id}: gap {ft.min() - top:.0f} Hz at {top:.0f} Hz -> stop")
        break
    specs.append((ft, asd))
    top = ft.max()
used = len(specs)
print(f"concatenated {used} runs, coverage {top/1e6:.3f} MHz")

def assemble(mode):
    """splice tiles; overlap kept by the low-f tile ('low') or overwritten by the high-f tile ('high')"""
    ftc, asdc = specs[0]
    for ft, asd in specs[1:]:
        if mode == "low":
            m = ft > ftc.max()
            ftc, asdc = np.concatenate([ftc, ft[m]]), np.concatenate([asdc, asd[m]])
        else:
            keep = ftc < ft.min()
            ftc, asdc = np.concatenate([ftc[keep], ft]), np.concatenate([asdc[keep], asd])
    return ftc, asdc

plt.rcParams.update({
    "font.size": 8, "axes.linewidth": 0.6, "legend.frameon": False,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})

def make_fig(x, y, label, fname, xscale="log", xlabel="frequency (Hz)"):
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax.plot(x, y, color="#1baf7a", lw=0.3, alpha=0.9, label=label)
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

# reference / filter-shape check
fig, ax = plt.subplots(figsize=(3.4, 2.6))
off = fref - 0.5 * (fref.min() + fref.max())
ax.plot(off, sref, lw=0.3, alpha=0.5, color="#898781", label=f"run {REF_RUN} raw")
ax.plot(off, ref, lw=1.0, color="#1baf7a", label=f"rolling avg (WIN={WIN})")
ax.set_yscale("log")
ax.set_xlabel("offset from window centre (Hz)"); ax.set_ylabel("V (arb.)")
ax.legend(fontsize=7, loc="lower center", handletextpad=0.4)
plt.tight_layout(pad=0.3)
plt.savefig(BASE + "_refshape.png", dpi=300, bbox_inches="tight")

# overlay of all compensated tiles (unspliced) -> overlap regions visible
fig, ax = plt.subplots(figsize=(3.4, 2.6))
for ft, asd in specs:
    ax.plot(ft, asd, lw=0.3, alpha=0.7)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("frequency (Hz)")
ax.set_ylabel(r"$\sqrt{S_I}$ (fA/$\sqrt{\mathrm{Hz}}$)")
ax.grid(True, which="major", color="#e1e0d9", lw=0.4)
ax.set_axisbelow(True)
ax.set_title(f"tile overlay ({used} runs, unspliced)", fontsize=7)
plt.tight_layout(pad=0.3)
plt.savefig(BASE + "_overlay.png", dpi=300, bbox_inches="tight")

for mode, suffix in [("low", ""), ("high", "_highwins")]:
    ftc, asdc = assemble(mode)
    if mode == "low":                          # remove the probe tone
        sel = np.abs(ftc - PROBE_F) < 500
        i = np.flatnonzero(sel)[np.argmax(asdc[sel])]
        keep = np.ones(asdc.size, bool)
        keep[i - GLITCH_HW:i + GLITCH_HW + 1] = False
        print(f"probe tone removed: {asdc[i]:.3g} fA/rtHz at {ftc[i]:.1f} Hz")
        ftc, asdc = ftc[keep], asdc[keep]
    lab = f"v47 baseband ({used} runs, {mode}-f tile wins)"
    make_fig(ftc, asdc, lab, BASE + suffix)
    make_fig(ftc / 1e3, asdc, lab, BASE + suffix + "_linx",
             xscale="linear", xlabel="frequency (kHz)")
plt.show()
