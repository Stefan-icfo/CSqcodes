#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ===================== USER CONFIG =====================
RUN_ID  = 346
DB_PATH = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v8.db"
OUT_DIR = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Figures"
os.makedirs(OUT_DIR, exist_ok=True)

# Integration windows (MHz)
PEAK_WINDOWS = {
    "Peak @142.099": (142.097, 142.101),
    "Peak @142.124": (142.122, 142.126),
}

# Tell the code which window is I^(+) and which is I^(-)
PLUS_WINDOW_NAME  = "Peak @142.124"   # choose one key from PEAK_WINDOWS
MINUS_WINDOW_NAME = "Peak @142.099"   # the other key from PEAK_WINDOWS
# =======================================================

qc.config["core"]["db_location"] = DB_PATH

# ---- helpers to load 1D PSD ----
CANDIDATE_SIGNAL_KEYS = [
    "avg_avg_psd_nodrive", "avg_psd_nodrive",
    "avg_avg_psd", "avg_psd",
    "V_fft_avg_avg", "PSD", "psd"
]
FREQ_KEY_CANDIDATES = ["freq_param", "frequency", "freq", "f"]

def load_1d_sorted(run_id):
    ds = qc.load_by_id(run_id)
    pd = ds.get_parameter_data()
    dep = next((k for k in CANDIDATE_SIGNAL_KEYS if k in pd), next(iter(pd.keys())))
    blk = pd[dep]
    sp_key = next(
        (k for k in blk if k != dep and any(t in k.lower() for t in FREQ_KEY_CANDIDATES)),
        next((k for k in blk if k != dep), None)
    )
    if sp_key is None:
        raise KeyError("No frequency setpoint key found.")

    x = np.asarray(blk[sp_key]).ravel() / 1e6  # Hz -> MHz
    y = np.asarray(blk[dep]).ravel()
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    idx = np.argsort(x)
    return x[idx], y[idx]

# ---- load data ----
x, y = load_1d_sorted(RUN_ID)
print(f"[Run {RUN_ID}] Loaded: N={len(x)}, f∈[{x.min():.6f},{x.max():.6f}] MHz")

# ---- integrate areas ----
def area_in_window(x, y, a, b):
    sel = (x >= a) & (x <= b)
    if np.count_nonzero(sel) < 2:
        raise RuntimeError(f"Not enough points in [{a},{b}] MHz.")
    return np.trapz(y[sel], x[sel])  # W/Hz * MHz

areas = {}
for name, (a, b) in PEAK_WINDOWS.items():
    A = area_in_window(x, y, a, b)
    areas[name] = A
    print(f"{name}: Area = {A:.3e} W/Hz·MHz  |  {A/1e-12:.3e} pW/Hz·MHz")

# ---- assign plus/minus and compute ratios ----
I_plus  = abs(areas[PLUS_WINDOW_NAME])
I_minus = abs(areas[MINUS_WINDOW_NAME])
ratio   = I_plus / I_minus if I_minus != 0 else np.inf

# Power ratio in dB (small/big and plus/minus both shown)
A_small, A_big = min(I_plus, I_minus), max(I_plus, I_minus)
ratio_small_big = A_small / A_big
squeezing_dB = 10 * np.log10(ratio_small_big)

print("\n=== Ratios ===")
print(f"I(+) = {I_plus:.3e}  |  I(-) = {I_minus:.3e}   (W/Hz·MHz)")
print(f"I(+)/I(-) = {ratio:.6f}")
print(f"Squeezing (power) = 10*log10(min/max) = {squeezing_dB:.2f} dB")

# ---- Huber et al. squeezing parameter φ from sideband asymmetry ----
# hi-branch:  I(+)/I(-) = 1 / tanh^2(phi_hi)  -> phi_hi = atanh( sqrt(1/ratio) )
# lo-branch:  I(+)/I(-) = tanh^2(phi_lo)      -> phi_lo = atanh( sqrt(ratio) )

def safe_atanh(z):
    z = np.clip(z, -0.999999999, 0.999999999)
    return 0.5 * np.log((1+z)/(1-z))

phi_hi = safe_atanh(np.sqrt(1.0/ratio)) if np.isfinite(ratio) and ratio>0 else np.nan
phi_lo = safe_atanh(np.sqrt(ratio))      if np.isfinite(ratio) and ratio>0 else np.nan

print("\n=== Squeezing parameter φ from Huber et al. ===")
print(f"Given PLUS='{PLUS_WINDOW_NAME}', MINUS='{MINUS_WINDOW_NAME}':")
print(f"  φ_hi  (I(+)/I(-)=1/tanh^2 φ) = {phi_hi:.6f}  rad = {np.degrees(phi_hi):.3f}°")
print(f"  φ_lo  (I(+)/I(-)=tanh^2  φ) = {phi_lo:.6f}  rad = {np.degrees(phi_lo):.3f}°")

# ---- plot with windows and annotations ----
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(x, y/1e-12, lw=1.2, color="tab:blue", label=f"Run {RUN_ID} PSD")

for name, (a, b) in PEAK_WINDOWS.items():
    sel = (x >= a) & (x <= b)
    lab = f"{name} area"
    ax.fill_between(x[sel], 0, y[sel]/1e-12, alpha=0.35, label=lab)

ax.set_xlabel("Frequency (MHz)", fontsize=16, fontname="Calibri")
ax.set_ylabel("PSD (pW/Hz)",   fontsize=16, fontname="Calibri")
ax.set_title(f"Run {RUN_ID}: areas, ratios and squeezing", fontsize=16, fontname="Calibri")
ax.grid(True, alpha=0.3)
ax.legend()

txt = (f"I(+)/I(-) = {ratio:.3e}\n"
       f"Squeezing = 10·log10(min/max) = {squeezing_dB:.2f} dB\n"
       f"φ_hi = {phi_hi:.3e} rad ({np.degrees(phi_hi):.2f}°)\n"
       f"φ_lo = {phi_lo:.3e} rad ({np.degrees(phi_lo):.2f}°)")
ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top",
        fontsize=12, bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

plt.tight_layout()
out_png = os.path.join(OUT_DIR, f"run_{RUN_ID}_areas_squeezing_huber.png")
plt.savefig(out_png, dpi=220)
print("Saved figure to:", out_png)
plt.show()


