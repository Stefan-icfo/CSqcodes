import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
from database import *

# paper figure: concatenated baseband spectra of v47 (Dec) AND v18 (Oct) on one linear axis.
# Same pipeline as paper_noisefloor_v47_concat; coverage holes become line breaks (NaN).
# NOTE: v18 uses the same GAIN/Z calibration as v47 (chain assumed unchanged Oct->Dec;
# supported by identical raw floors).

F_MIX, RBW, WIN, GLITCH_HW = 2.5e6, 1.676, 1001, 1
FMIN = 1.0
GAIN = 200 * 5.65   # amplification chain
Z    = 7.5e3        # Ohm, V -> I conversion
PROBE_F = 100e3     # Hz, probe tone: removed if a clear spike sits there

DATASETS = [
    dict(label="v47 (8 Dec)", color="#1baf7a",
         db=r'D:\databases CD12_B5_F4\CD12_B5_F4v47_08_12_25.db', ref=72,
         runs=[78, 80, 82, 84, 86, 88, 90, 93, 95, 97, 99, 101, 103, 105, 107, 109,
               111, 114, 116, 118, 9]),                     # contiguous 0..1.055 MHz
    dict(label="v18 (18 Oct)", color="#eda100",
         db=r'D:\databases CD12_B5_F4\CD12_B5_F4v18_171025.db', ref=198,
         runs=[178, 180, 182]),                             # true 0..155 kHz + 0.945..1.055 MHz
]

def get_spec(run_id):
    f, s = extract_1d(run_id, data_1d_name="V_fft_avg_avg",
                      setpoint_name="freq_param", plot=False, return_exp_name=False)
    return np.asarray(f, float), np.asarray(s, float)

def rolling_mean(y, win):
    pad = win // 2
    return np.convolve(np.pad(y, pad, mode="edge"), np.ones(win) / win, mode="valid")

def build(ds):
    qc.config["core"]["db_location"] = ds["db"]
    fref, sref = get_spec(ds["ref"])
    clean = sref.copy()
    i0 = np.argmax(clean)
    clean[max(0, i0 - GLITCH_HW):i0 + GLITCH_HW + 1] = np.nan
    bad = np.isnan(clean)
    clean[bad] = np.interp(fref[bad], fref[~bad], clean[~bad])
    ref_norm = rolling_mean(clean, WIN)
    ref_norm /= ref_norm.max()

    ftc, asdc = np.array([]), np.array([])
    for run_id in ds["runs"]:
        f, s = get_spec(run_id)
        R = ref_norm if s.size == ref_norm.size else np.interp(
            np.linspace(0, 1, s.size), np.linspace(0, 1, ref_norm.size), ref_norm)
        ft = f - F_MIX
        m = ft >= FMIN
        ft, asd = ft[m], (s[m] / R[m]) / np.sqrt(RBW) / (GAIN * Z) * 1e15   # fA/sqrt(Hz)
        if ftc.size and ft.min() > ftc[-1]:                 # coverage hole -> break the line
            ftc = np.append(ftc, 0.5 * (ftc[-1] + ft.min()))
            asdc = np.append(asdc, np.nan)
        m2 = ft > (ftc[np.isfinite(asdc)][-1] if ftc.size else -np.inf)  # low-f tile wins
        ftc, asdc = np.concatenate([ftc, ft[m2]]), np.concatenate([asdc, asd[m2]])
    # remove the probe tone if a clear spike sits at PROBE_F
    sel = (np.abs(ftc - PROBE_F) < 500) & np.isfinite(asdc)
    if sel.any():
        i = np.flatnonzero(sel)[np.nanargmax(asdc[sel])]
        if asdc[i] > 5 * np.nanmedian(asdc[sel]):
            keep = np.ones(asdc.size, bool)
            keep[i - GLITCH_HW:i + GLITCH_HW + 1] = False
            print(f"{ds['label']}: probe tone removed, {asdc[i]:.3g} fA/rtHz at {ftc[i]:.1f} Hz")
            ftc, asdc = ftc[keep], asdc[keep]
    print(f"{ds['label']}: {len(ds['runs'])} runs, coverage to {np.nanmax(ftc)/1e6:.3f} MHz")
    return ftc, asdc

plt.rcParams.update({
    "font.size": 8, "axes.linewidth": 0.6, "legend.frameon": False,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
import os
FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)
BASE = os.path.join(FIGDIR, os.path.basename(__file__).replace(".py", ""))

fig, ax = plt.subplots(figsize=(3.4, 2.6))
for ds in DATASETS:
    ftc, asdc = build(ds)
    ax.plot(ftc / 1e3, asdc, lw=0.3, color=ds["color"], alpha=0.85, label=ds["label"])
ax.set_yscale("log")
ax.set_xlabel("frequency (kHz)")
ax.set_ylabel(r"$\sqrt{S_I}$ (fA/$\sqrt{\mathrm{Hz}}$)")
ax.grid(True, which="major", color="#e1e0d9", lw=0.4)
ax.set_axisbelow(True)
ax.legend(fontsize=7, loc="upper right", handletextpad=0.4)
plt.tight_layout(pad=0.3)
plt.savefig(BASE + ".png", dpi=300, bbox_inches="tight")
plt.savefig(BASE + ".pdf", bbox_inches="tight")
plt.show()
