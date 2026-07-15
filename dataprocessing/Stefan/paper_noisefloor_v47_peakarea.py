import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
from database import *

# based on paper_noisefloor_v47_concat: subtract the flat background level from the
# composite PSD; everything above it = peaks (charge noise). Compare integrated powers.

F_MIX, RBW, WIN, GLITCH_HW = 2.5e6, 1.676, 1001, 1
FMIN = 1.0
GAIN = 200 * 5.65   # amplification chain
Z    = 7.5e3        # Ohm, V -> I conversion
PROBE_F = 100e3     # Hz, probe tone (not charge noise) -> removed
DB, REF_RUN = r'D:\databases CD12_B5_F4\CD12_B5_F4v47_08_12_25.db', 72
CANDIDATES = [78, 80, 82, 84, 86, 88, 90, 93, 95, 97, 99, 101, 103, 105, 107, 109,
              111, 114, 116, 118, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27]

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
ref_norm = rolling_mean(clean, WIN)
ref_norm /= ref_norm.max()

def spec_comp(run_id):
    f, s = get_spec(run_id)
    R = ref_norm if s.size == ref_norm.size else np.interp(
        np.linspace(0, 1, s.size), np.linspace(0, 1, ref_norm.size), ref_norm)
    ft = f - F_MIX
    m = ft >= FMIN
    return ft[m], (s[m] / R[m]) / np.sqrt(RBW) / (GAIN * Z)   # A/sqrt(Hz)

specs, top = [], None
for run_id in CANDIDATES:
    ft, asd = spec_comp(run_id)
    if top is not None and ft.min() > top:
        break
    specs.append((ft, asd))
    top = ft.max()

ftc, asdc = specs[0]
for ft, asd in specs[1:]:                       # low-f tile wins in overlaps
    m = ft > ftc.max()
    ftc, asdc = np.concatenate([ftc, ft[m]]), np.concatenate([asdc, asd[m]])
sel = np.abs(ftc - PROBE_F) < 500               # remove the probe tone
i = np.flatnonzero(sel)[np.argmax(asdc[sel])]
keep = np.ones(asdc.size, bool)
keep[i - GLITCH_HW:i + GLITCH_HW + 1] = False
ftc, asdc = ftc[keep], asdc[keep]

# ---- flat background level; everything above it = peaks ----
S  = asdc ** 2                                  # PSD, A^2/Hz
bg = np.median(S)                               # flat background level
df = np.gradient(ftc)
P_bg    = bg * np.sum(df)                       # A^2
P_peaks = np.sum((S - bg) * df)                 # all the rest
P_tot   = np.sum(S * df)
print(f"band {ftc.min():.0f} Hz .. {ftc.max()/1e6:.3f} MHz")
print(f"background level: {np.sqrt(bg)*1e15:.1f} fA/rtHz")
print(f"background power: {P_bg:.3e} A^2  (rms {np.sqrt(P_bg)*1e12:.1f} pA)")
print(f"peak power      : {P_peaks:.3e} A^2  (rms {np.sqrt(P_peaks)*1e12:.1f} pA)")
print(f"peaks/background = {P_peaks/P_bg:.2f}   (total {P_tot:.3e} A^2)")

plt.rcParams.update({
    "font.size": 8, "axes.linewidth": 0.6, "legend.frameon": False,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
import os
FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)
BASE = os.path.join(FIGDIR, os.path.basename(__file__).replace(".py", ""))

# fig 1: PSD with the flat background level
fig, ax = plt.subplots(figsize=(3.4, 2.6))
ax.plot(ftc / 1e3, S, lw=0.3, color="#1baf7a", alpha=0.8, label="composite PSD")
ax.axhline(bg, lw=0.8, color="#2a78d6", label="flat background level")
ax.set_yscale("log")
ax.set_xlabel("frequency (kHz)")
ax.set_ylabel(r"$S_I$ (A$^2$/Hz)")
ax.grid(True, which="major", color="#e1e0d9", lw=0.4)
ax.set_axisbelow(True)
ax.legend(fontsize=6, loc="upper right", handletextpad=0.4)
plt.tight_layout(pad=0.3)
plt.savefig(BASE + "_decomp.png", dpi=300, bbox_inches="tight")
plt.savefig(BASE + "_decomp.pdf", bbox_inches="tight")

# fig 2: cumulative integrated power of the two components
fig, ax = plt.subplots(figsize=(3.4, 2.6))
cum_bg = np.cumsum(bg * df) * 1e24              # pA^2
cum_pk = np.cumsum((S - bg) * df) * 1e24
ax.plot(ftc / 1e3, cum_bg, lw=1.0, color="#2a78d6", label="background")
ax.plot(ftc / 1e3, cum_pk, lw=1.0, color="#e34948", label="peaks (rest)")
ax.set_xlabel("frequency (kHz)")
ax.set_ylabel(r"cumulative power (pA$^2$)")
ax.grid(True, which="major", color="#e1e0d9", lw=0.4)
ax.set_axisbelow(True)
ax.legend(fontsize=7, loc="upper left", handletextpad=0.4)
plt.tight_layout(pad=0.3)
plt.savefig(BASE + "_cumulative.png", dpi=300, bbox_inches="tight")
plt.savefig(BASE + "_cumulative.pdf", bbox_inches="tight")
plt.show()
