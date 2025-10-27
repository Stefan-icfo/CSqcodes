#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Steps 1–4: Mass-free displacement calibration (thermal + driven), PSD kept in W/Hz (50 Ohm system).

Step 1: Background subtract thermal PSD (W/Hz) and integrate -> <P> [W], then <V^2> = R * <P> [V^2].
Step 2: Estimate Q from the thermal peak FWHM (uses background-subtracted PSD).
Step 3: Read driven sweep and get V_dr at resonance (near f0).
Step 4: Compute mass-free gain g = Q * F_w * <V^2> / (k_B * T * V_dr(omega0))  [V/m].

Notes:
- NO unit conversion is applied to arrays: PSD stays in W/Hz.
- R = 50 Ohm is used only to report <V^2> from the integrated power.
- For the drive chain, we convert the instrument voltage to the device via the given attenuation (dB).

This version avoids non-ASCII symbols in print strings/titles to prevent UnicodeEncodeError on Windows consoles.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ======================= USER CONFIG =======================
DB_PATH = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v11.db"

# Data runs
THERMAL_RUN_ID    = 1623   # thermal spectrum (data)
BACKGROUND_RUN_ID = 1606   # background spectrum
DRIVEN_RUN_ID     = 1509   # driven sweep (voltage vs frequency); set None to skip

# Subtraction
SUBTRACTION_MODE   = "interpolate"  # "index" or "interpolate" (recommended)
VOLTAGE_DOMAIN_SUB = False          # keep pure power subtraction: P_corr = P_data - P_bg  (still W/Hz)

# Physics constants / environment
R_OHM         = 50.0         # RF chain impedance for reporting <V^2> = R * <P>
K_B           = 1.380649e-23 # Boltzmann [J/K]
TEMPERATURE_K = 0.035        # <-- set your mode temperature [K]

# Drive chain and force model (electrostatic small-signal: F ≈ (dC/dx) * Vdc * v_ac )
V_SRC_VPK      = 37.5e-3     # instrument output amplitude [V_peak] at the source (given)
ATTEN_DB       = 43.0        # total attenuation from source to device [dB] (given)
EXTRA_LINE_GAIN = 1.0        # multiply by any additional known line factor if needed (typically 1)
V_DC_GATE      = 0.5         # <-- set your DC gate voltage [V]
DCDX_F_PER_M   = 1.0e-12     # <-- set your dC/dx [F/m], example value! put your number

# Driven-peak selection
DRIVEN_WINDOW_MULT = 2.0     # search ±(DRIVEN_WINDOW_MULT * Gamma) around f0 for V_dr(omega0)

# Plot / output
USE_INDEX_AXIS_FOR_PLOT = False
SHOW_FIGS = True
OUT_DIR = os.path.abspath("./calib_mass_free_W_per_Hz")
os.makedirs(OUT_DIR, exist_ok=True)
print("Output folder:", OUT_DIR)
try:
    os.startfile(OUT_DIR)  # Windows convenience
except Exception:
    pass

# (Optional) force specific keys if auto-detect picks the wrong one
FORCE_DEP_KEY  = None          # e.g. "avg_avg_psd_nodrive"
FORCE_FREQ_KEY = None

# Matplotlib font (safe)
mpl.rcParams["font.family"] = "Calibri"

# Likely parameter names
CANDIDATE_SIGNAL_KEYS = [
    "avg_avg_psd_nodrive", "avg_psd_nodrive",
    "avg_avg_psd", "avg_psd",
    "PSD", "psd", "Svv", "V_fft_avg_avg", "V_r", "Amplitude"
]
FREQ_KEY_CANDIDATES = ["freq_param", "compressed_freq_real", "freq", "frequency", "Frequency", "f"]

# ---------------------- Helpers ----------------------

def save_csv(path, header, *cols):
    arr = np.column_stack(cols)
    np.savetxt(path, arr, delimiter=",", header=header, comments="")
    print("[saved]", path)

def load_trace_1d(run_id):
    """Load one dependent 1D trace (keeps acquisition order) and, if present, its frequency setpoint."""
    ds = load_by_id(run_id)
    pd = ds.get_parameter_data()

    if FORCE_DEP_KEY is not None:
        dep = FORCE_DEP_KEY
        if dep not in pd:
            raise KeyError(f"Forced dep key '{dep}' not found in dataset.")
    else:
        dep = next((k for k in CANDIDATE_SIGNAL_KEYS if k in pd), next(iter(pd.keys())))
    blk = pd[dep]

    if FORCE_FREQ_KEY is not None:
        fkey = FORCE_FREQ_KEY
        if fkey not in blk:
            raise KeyError(f"Forced freq key '{fkey}' not found in block for dep '{dep}'.")
    else:
        fkey = None
        for k in blk.keys():
            if k == dep:
                continue
            if any(tag.lower() in k.lower() for tag in FREQ_KEY_CANDIDATES):
                fkey = k
                break

    y = np.asarray(blk[dep]).ravel().astype(float)
    f = np.asarray(blk[fkey]).ravel().astype(float) if fkey is not None else None

    m = np.isfinite(y)
    if f is not None:
        m &= np.isfinite(f)
        f = f[m]
    y = y[m]

    print(f"[run {run_id}] dep='{dep}', fkey='{fkey}', N={len(y)}")
    return f, y, dep, fkey

def subtract_background(thermal_id, bg_id, mode="interpolate", voltage_domain=False):
    """Return f_Hz, P_data(W/Hz), P_bg_aligned(W/Hz), P_corr(W/Hz)."""
    f_d, y_d, _, _ = load_trace_1d(thermal_id)
    f_b, y_b, _, _ = load_trace_1d(bg_id)

    P_d = np.asarray(y_d, dtype=float)  # W/Hz
    P_b = np.asarray(y_b, dtype=float)  # W/Hz

    if mode.lower().startswith("index"):
        n = min(len(P_d), len(P_b))
        Pd = P_d[:n]; Pb = P_b[:n]
        if voltage_domain:
            Pc = (np.sqrt(np.clip(Pd,0,None)) - np.sqrt(np.clip(Pb,0,None)))**2
        else:
            Pc = Pd - Pb
        f_used = (None if (f_d is None or USE_INDEX_AXIS_FOR_PLOT) else f_d[:n])
        return f_used, Pd[:n], Pb[:n], Pc

    if f_d is None or f_b is None:
        raise ValueError("Mode 'interpolate' requires both frequency axes.")
    Pbg_interp = np.interp(f_d, f_b, P_b, left=P_b[0], right=P_b[-1])
    if voltage_domain:
        Pc = (np.sqrt(np.clip(P_d,0,None)) - np.sqrt(np.clip(Pbg_interp,0,None)))**2
    else:
        Pc = P_d - Pbg_interp
    return f_d, P_d, Pbg_interp, Pc

def integrate_psd_W_per_Hz(Pcorr, f_Hz, one_sided=True, clip_neg=True):
    """Return <P> [W] from S_PP(f) [W/Hz]."""
    Pc = np.asarray(Pcorr, dtype=float)
    f  = np.asarray(f_Hz, dtype=float)
    if clip_neg:
        Pc = np.clip(Pc, 0.0, None)
    m = np.isfinite(Pc) & np.isfinite(f)
    if np.count_nonzero(m) < 2:
        return np.nan
    P_total = np.trapz(Pc[m], f[m])    # W
    if not one_sided:
        P_total *= 2.0
    return P_total

def estimate_f0_Q_from_halfmax(f_Hz, Pcorr_W_per_Hz):
    """Estimate f0 and Q from the background-subtracted thermal peak (assumed Lorentzian)."""
    f = np.asarray(f_Hz, dtype=float)
    S = np.asarray(Pcorr_W_per_Hz, dtype=float)
    m = np.isfinite(f) & np.isfinite(S)
    f = f[m]; S = S[m]
    if len(S) < 5:
        return np.nan, np.nan, np.nan

    i0 = int(np.nanargmax(S))
    f0 = float(f[i0])
    S0 = float(S[i0])
    if S0 <= 0:
        return np.nan, np.nan, np.nan
    half = S0 / 2.0

    # find left crossing
    iL = i0
    while iL > 0 and S[iL] > half:
        iL -= 1
    # find right crossing
    iR = i0
    nmax = len(S) - 1
    while iR < nmax and S[iR] > half:
        iR += 1

    if iL == iR:
        return f0, np.nan, np.nan

    # linear interpolate for better FWHM
    def interp_x(i1, i2):
        x1, x2 = f[i1], f[i2]
        y1, y2 = S[i1], S[i2]
        if y2 == y1:
            return x1
        return x1 + (half - y1) * (x2 - x1) / (y2 - y1)

    fL = interp_x(iL, min(iL+1, len(S)-1))
    fR = interp_x(max(iR-1, 0), iR)
    FWHM = float(abs(fR - fL))  # Gamma in Hz
    Q = f0 / FWHM if FWHM > 0 else np.nan
    return f0, FWHM, Q

def db_to_linear_voltage(att_db):
    """Convert dB attenuation to linear voltage ratio."""
    return 10.0 ** (-att_db / 20.0)

def compute_drive_force(V_src_vpk, atten_db, extra_gain, Vdc, dCdx_F_per_m):
    """
    Electrostatic small-signal force: F_w ≈ (dC/dx) * Vdc * v_ac ,
    with v_ac = V_src_vpk * 10^(-atten/20) * extra_gain
    Returns F_w in newtons and the device voltage v_ac [V_peak].
    """
    v_ac_device = V_src_vpk * db_to_linear_voltage(atten_db) * extra_gain
    Fw = dCdx_F_per_m * Vdc * v_ac_device
    return Fw, v_ac_device

def get_Vdr_at_resonance(driven_run_id, f0_Hz, Gamma_Hz, window_mult=2.0):
    """
    From the driven sweep, get V_dr(omega0): pick the maximum absolute voltage within ±(window_mult * Gamma) around f0.
    Returns (V_dr_at_res, f_at_Vdr).
    """
    f, y, dep, fkey = load_trace_1d(driven_run_id)
    if f is None or len(f) == 0:
        raise RuntimeError("Driven run has no frequency axis.")
    V = np.asarray(y, dtype=float)

    # window selection
    fmin = f0_Hz - window_mult * Gamma_Hz
    fmax = f0_Hz + window_mult * Gamma_Hz
    mask = (f >= fmin) & (f <= fmax)
    if np.count_nonzero(mask) < 3:
        # fallback: nearest point to f0
        idx = int(np.argmin(np.abs(f - f0_Hz)))
        return float(np.abs(V[idx])), float(f[idx])

    idx_local = np.argmax(np.abs(V[mask]))
    Vdr = float(np.abs(V[mask][idx_local]))
    f_at = float(f[mask][idx_local])
    return Vdr, f_at

def plot_thermal(f, Pdat, Pbg, Pcorr, tag):
    # corrected
    plt.figure(figsize=(10, 6))
    x = (f / 1e6) if f is not None else np.arange(len(Pcorr))
    xlabel = "Frequency (MHz)" if f is not None else "Point index"
    plt.plot(x, Pcorr, "-", lw=1.8, label="Thermal (BG-subtracted)")
    plt.xlabel(xlabel); plt.ylabel("PSD [W/Hz]"); plt.title(f"Thermal PSD after BG subtraction {tag}")
    plt.legend(); out = os.path.join(OUT_DIR, f"thermal_corrected{tag.replace(' ', '_')}.png")
    plt.tight_layout(); plt.savefig(out, dpi=200); print("[saved]", out)
    if SHOW_FIGS: plt.show(); plt.close()

    # diagnostic
    plt.figure(figsize=(10, 6))
    if f is None:
        n = len(Pcorr); x = np.arange(n); xlabel = "Point index"
        plt.plot(x, Pdat[:n], alpha=0.7, label="Data")
        plt.plot(x, Pbg[:n],  alpha=0.7, label="Background (aligned)")
        plt.plot(x, Pcorr, lw=1.8, label="Corrected")
    else:
        x = f / 1e6; xlabel = "Frequency (MHz)"
        plt.plot(x, Pdat, alpha=0.7, label="Data")
        plt.plot(x, Pbg,  alpha=0.7, label="Background (aligned)")
        plt.plot(x, Pcorr, lw=1.8, label="Corrected")
    plt.xlabel(xlabel); plt.ylabel("PSD [W/Hz]"); plt.title("Thermal: data vs background vs corrected")
    plt.legend(); out = os.path.join(OUT_DIR, f"thermal_data_bg_corrected_{THERMAL_RUN_ID}_{BACKGROUND_RUN_ID}.png")
    plt.tight_layout(); plt.savefig(out, dpi=200); print("[saved]", out)
    if SHOW_FIGS: plt.show(); plt.close()

# ------------------------ Main ------------------------

def main():
    initialise_or_create_database_at(DB_PATH)

    # ---- Step 1: thermal background subtraction & integral ----
    f, Pdat, Pbg, Pcorr = subtract_background(
        THERMAL_RUN_ID, BACKGROUND_RUN_ID,
        mode=SUBTRACTION_MODE,
        voltage_domain=VOLTAGE_DOMAIN_SUB
    )
    if f is None:
        raise RuntimeError("Thermal run has no frequency axis; cannot integrate in physical units.")

    tag = f"(run {THERMAL_RUN_ID} - {BACKGROUND_RUN_ID})"
    plot_thermal(f, Pdat, Pbg, Pcorr, tag)

    P_total_W = integrate_psd_W_per_Hz(Pcorr, f, one_sided=True, clip_neg=True)
    V2_total  = R_OHM * P_total_W
    print("[Step 1] <P>   = {:.6e} W".format(P_total_W))
    print("[Step 1] <V^2> = {:.6e} V^2   (R = {:.0f} Ohm)".format(V2_total, R_OHM))

    # ---- Step 2: estimate f0, Gamma, Q from the thermal peak ----
    f0_Hz, Gamma_Hz, Q_est = estimate_f0_Q_from_halfmax(f, Pcorr)
    print("[Step 2] f0 = {:.6f} Hz,  Gamma(FWHM) = {:.6f} Hz,  Q = {:.3g}".format(f0_Hz, Gamma_Hz, Q_est))

    # ---- Step 3: get driven voltage at resonance ----
    if DRIVEN_RUN_ID is None:
        raise RuntimeError("DRIVEN_RUN_ID is None; cannot complete Step 3/4.")
    Vdr_at_res_V, f_at_Vdr = get_Vdr_at_resonance(DRIVEN_RUN_ID, f0_Hz, Gamma_Hz, window_mult=DRIVEN_WINDOW_MULT)
    print("[Step 3] V_dr(omega0) ~= {:.6e} V at f = {:.6f} Hz".format(Vdr_at_res_V, f_at_Vdr))

    # ---- Drive force amplitude from the chain ----
    Fw_N, v_ac_device_V = compute_drive_force(V_SRC_VPK, ATTEN_DB, EXTRA_LINE_GAIN, V_DC_GATE, DCDX_F_PER_M)
    print("[Step 3] v_ac at device = {:.6e} V_peak  (source {:.3e} Vpk, atten {:.1f} dB)".format(
        v_ac_device_V, V_SRC_VPK, ATTEN_DB))
    print("[Step 3] F_w (model)    = {:.6e} N  using F ~= (dC/dx)*Vdc*v_ac".format(Fw_N))

    # ---- Step 4: mass-free calibration ----
    g_V_per_m = (Q_est * Fw_N * V2_total) / (K_B * TEMPERATURE_K * Vdr_at_res_V)
    print("[Step 4] g = {:.6e} V/m   (mass-free)".format(g_V_per_m))

    # ---- Save a small report and CSVs ----
    save_csv(os.path.join(OUT_DIR, "thermal_corrected_{}_minus_{}.csv".format(THERMAL_RUN_ID, BACKGROUND_RUN_ID)),
             "f_Hz,Pcorr_W_per_Hz", f, Pcorr)
    save_csv(os.path.join(OUT_DIR, "mass_free_report.csv"),
             "f0_Hz,Gamma_Hz,Q,<P>_W,<V2>_V2,Vdr_V,Fw_N,T_K,g_V_per_m",
             np.array([f0_Hz]), np.array([Gamma_Hz]), np.array([Q_est]),
             np.array([P_total_W]), np.array([V2_total]),
             np.array([Vdr_at_res_V]), np.array([Fw_N]),
             np.array([TEMPERATURE_K]), np.array([g_V_per_m])
    )

    # Human-readable TXT (UTF-8 so you can use symbols later if desired)
    report_txt = os.path.join(OUT_DIR, "mass_free_report.txt")
    with open(report_txt, "w", encoding="utf-8") as fh:
        fh.write("Mass-free calibration report\n")
        fh.write("DB: {}\n".format(DB_PATH))
        fh.write("Thermal run: {} | Background run: {} | Driven run: {}\n".format(
            THERMAL_RUN_ID, BACKGROUND_RUN_ID, DRIVEN_RUN_ID))
        fh.write("R = {:.1f} Ohm,  T = {:.6f} K\n".format(R_OHM, TEMPERATURE_K))
        fh.write("Drive: V_src = {:.6e} Vpk, atten = {:.1f} dB, Vdc = {:.6e} V, dC/dx = {:.6e} F/m\n".format(
            V_SRC_VPK, ATTEN_DB, V_DC_GATE, DCDX_F_PER_M))
        fh.write("v_ac(device) = {:.6e} Vpk,  F_w = {:.6e} N\n".format(v_ac_device_V, Fw_N))
        fh.write("Step 1: <P> = {:.6e} W,  <V^2> = {:.6e} V^2\n".format(P_total_W, V2_total))
        fh.write("Step 2: f0 = {:.6f} Hz, Gamma(FWHM) = {:.6f} Hz, Q = {:.6g}\n".format(f0_Hz, Gamma_Hz, Q_est))
        fh.write("Step 3: V_dr(omega0) = {:.6e} V at f = {:.6f} Hz\n".format(Vdr_at_res_V, f_at_Vdr))
        fh.write("Step 4: g = {:.6e} V/m (mass-free)\n".format(g_V_per_m))
    print("[saved]", report_txt)

    print("Done. Files in:", OUT_DIR)

if __name__ == "__main__":
    main()


