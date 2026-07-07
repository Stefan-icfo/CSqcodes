import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
from database import *

# paper figure: broadband noise-floor surveys v47 + v49 as point series.
# each point = MEDIAN compensated ASD of one spectrum (robust to sharp spikes);
# spectra whose window is wide relative to its centre get one median per log bin.

F_MIX   = 2.5e6      # Hz, mix-down tone; true freq = freq_param - F_MIX
RBW     = 1.676      # Hz
WIN     = 51         # reference smoothing window
GLITCH_HW = 1
FMIN    = 1.0        # Hz, discard bins below this
NBIN_DEC  = 4        # log bins per decade for wide (low-freq) windows
MIN_SAMP  = 8        # min samples per log bin
REL_W     = 0.5      # window/centre ratio above which a spectrum is "wide"
EDGE_TRIM = 0.4      # wide windows: keep |offset| <= EDGE_TRIM*width (edges compensate poorly)

V49 = dict(db=r'D:\databases CD12_B5_F4\CD12_B5_F4v49_15_12_25.db', ref=188,
           runs=[152, 155, 158, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 191, 194])
V47 = dict(db=r'D:\databases CD12_B5_F4\CD12_B5_F4v47_08_12_25.db', ref=72,
           runs=[78, 80, 82, 84, 86, 88, 90, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111,
                 114, 116, 118, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 30, 32, 34, 36, 38,
                 40, 42, 44, 46, 48, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 72])

def get_spec(run_id):
    f, s = extract_1d(run_id, data_1d_name="V_fft_avg_avg",
                      setpoint_name="freq_param", plot=False, return_exp_name=False)
    return np.asarray(f, float), np.asarray(s, float)

def rolling_mean(y, win):
    pad = win // 2
    return np.convolve(np.pad(y, pad, mode="edge"), np.ones(win) / win, mode="valid")

def build_ref(run_id):
    fref, sref = get_spec(run_id)
    clean = sref.copy()
    i0 = np.argmax(clean)
    clean[max(0, i0 - GLITCH_HW):i0 + GLITCH_HW + 1] = np.nan
    bad = np.isnan(clean)
    clean[bad] = np.interp(fref[bad], fref[~bad], clean[~bad])
    ref = rolling_mean(clean, WIN)
    return ref / ref.max()

def load_raw(cfg):
    """narrow spectra -> one median point each (bin-shift invariant);
    wide (low-freq) spectra -> pooled samples for log binning"""
    qc.config["core"]["db_location"] = cfg["db"]
    ref_norm = build_ref(cfg["ref"])
    narrow, pool_f, pool_a = [], [], []
    for run_id in cfg["runs"]:
        f, s = get_spec(run_id)
        R = ref_norm if s.size == ref_norm.size else np.interp(
            np.linspace(0, 1, s.size), np.linspace(0, 1, ref_norm.size), ref_norm)
        asd = (s / R) / np.sqrt(RBW)
        ft = f - F_MIX
        c, W = 0.5 * (ft.min() + ft.max()), ft.max() - ft.min()
        if c > 0 and W / c < REL_W:
            narrow.append((c, np.median(asd)))
        else:
            m = (ft >= FMIN) & (np.abs(ft - c) <= EDGE_TRIM * W)
            pool_f.append(ft[m]); pool_a.append(asd[m])
    return narrow, np.concatenate(pool_f), np.concatenate(pool_a)

def bin_pool(pool_f, pool_a, delta=0.0, phase=0.0):
    """one median per global log bin; grid shifted by delta [Hz] and/or phase [bins];
    median = robust to spikes"""
    fs = pool_f + delta
    m = fs >= FMIN
    k = np.floor(np.log10(fs[m]) * NBIN_DEC - phase).astype(int)
    pa = pool_a[m]
    out = {}
    for kb in np.unique(k):
        sel = k == kb
        if sel.sum() >= MIN_SAMP:
            out[kb] = (10 ** ((kb + phase + 0.5) / NBIN_DEC), np.median(pa[sel]))
    return out

def make_points(raw, delta=0.0, phase=0.0):
    pts = list(raw[0]) + list(bin_pool(raw[1], raw[2], delta, phase).values())
    pts.sort()
    return np.array(pts).T

raw47 = load_raw(V47)
raw49 = load_raw(V49)

GRID_PHASE = 0.33   # bin-grid offset (fraction of a bin), chosen via bin-shift robustness checks

plt.rcParams.update({
    "font.size": 8, "axes.linewidth": 0.6, "legend.frameon": False,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
pts49 = make_points(raw49, phase=GRID_PHASE)
pts47 = make_points(raw47, phase=GRID_PHASE)
fig, ax = plt.subplots(figsize=(3.4, 2.6))
ax.plot(*pts49, "o-", color="#2a78d6", ms=3.2, lw=0.8, alpha=0.9,
        markeredgewidth=0, label="v49 (26 Dec)")
ax.plot(*pts47, "s-", color="#1baf7a", ms=3.0, lw=0.8, alpha=0.9,
        markeredgewidth=0, label="v47 (8 Dec)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("frequency (Hz)")
ax.set_ylabel(r"ASD (V/$\sqrt{\mathrm{Hz}}$)")
ax.grid(True, which="major", color="#e1e0d9", lw=0.4)
ax.set_axisbelow(True)
ax.legend(fontsize=7, loc="upper right", handletextpad=0.4)
plt.tight_layout(pad=0.3)
import os
FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)
BASE = os.path.join(FIGDIR, os.path.basename(__file__).replace(".py", ""))
plt.savefig(BASE + ".png", dpi=300, bbox_inches="tight")
plt.savefig(BASE + ".pdf", bbox_inches="tight")
plt.show()
