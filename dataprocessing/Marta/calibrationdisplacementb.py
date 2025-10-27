#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Thermal background subtraction (keeping PSD in W/Hz) + integral:
  <P> = ∫ S_PP(f) df  [W], and <V^2> = R * <P> [V^2]
Driven sweep plotted as stored (no calibration).

No unit conversion applied to arrays; we only *report* <V^2> using R=50 Ω.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ======================= USER CONFIG =======================
DB_PATH = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v11.db"

THERMAL_RUN_ID    = 1623   # thermal spectrum (data)
BACKGROUND_RUN_ID = 1606   # background spectrum
DRIVEN_RUN_ID     = 1509   # driven sweep (voltage vs frequency); set None to skip

SUBTRACTION_MODE   = "interpolate"  # "index" or "interpolate" (recommended)
VOLTAGE_DOMAIN_SUB = False          # keep pure power subtraction: P_corr = P_data - P_bg  (still W/Hz)

USE_INDEX_AXIS_FOR_PLOT = False
SHOW_FIGS = True

OUT_DIR = os.path.abspath("./fig_subtract_integrate_W_per_Hz")
os.makedirs(OUT_DIR, exist_ok=True)
print("Output folder:", OUT_DIR)
try:
    os.startfile(OUT_DIR)  # Windows convenience
except Exception:
    pass

# (Optional) force specific keys if auto-detect picks the wrong one
FORCE_DEP_KEY  = None          # e.g. "avg_avg_psd_nodrive"
FORCE_FREQ_KEY = None

# ======= UNITS (display only; arrays remain as stored) =======
PSD_UNITS = "W/Hz"     # your PSD is power spectral density
R_OHM     = 50.0       # impedance used only to *report* <V^2> = R * <P>
# ============================================================

mpl.rcParams["font.family"] = "Calibri"

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
    return f_d, P_d, Pbg_interp, Pc  # all in W/Hz

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

def plot_thermal_corrected(f, P_corr, title_suffix=""):
    plt.figure(figsize=(10, 6))
    if f is None or USE_INDEX_AXIS_FOR_PLOT:
        x = np.arange(len(P_corr)); xlabel = "Point index"
    else:
        x = f / 1e6; xlabel = "Frequency (MHz)"
    plt.plot(x, P_corr, "-", lw=1.8, label="Thermal (BG-subtracted)")
    plt.xlabel(xlabel)
    plt.ylabel("PSD [W/Hz]")
    plt.title(f"Thermal PSD after BG subtraction {title_suffix}")
    plt.legend()
    out = os.path.join(OUT_DIR, f"thermal_corrected{title_suffix.replace(' ', '_')}.png")
    plt.tight_layout(); plt.savefig(out, dpi=200); print("[saved]", out)
    if SHOW_FIGS: plt.show()
    plt.close()

def plot_thermal_diagnostic(f, Pdat, Pbg, Pcorr, thermal_id, bg_id):
    plt.figure(figsize=(10, 6))
    if f is None or USE_INDEX_AXIS_FOR_PLOT:
        n = len(Pcorr); x = np.arange(n); xlabel = "Point index"
        plt.plot(x, Pdat[:n], alpha=0.7, label="Data")
        plt.plot(x, Pbg[:n],  alpha=0.7, label="Background (aligned)")
        plt.plot(x, Pcorr, lw=1.8, label="Corrected")
    else:
        x = f / 1e6; xlabel = "Frequency (MHz)"
        plt.plot(x, Pdat, alpha=0.7, label="Data")
        plt.plot(x, Pbg,  alpha=0.7, label="Background (aligned)")
        plt.plot(x, Pcorr, lw=1.8, label="Corrected")
    plt.xlabel(xlabel); plt.ylabel("PSD [W/Hz]")
    plt.title("Thermal PSD: data vs background vs corrected (W/Hz)")
    plt.legend()
    out = os.path.join(OUT_DIR, f"thermal_data_bg_corrected_{thermal_id}_{bg_id}.png")
    plt.tight_layout(); plt.savefig(out, dpi=200); print("[saved]", out)
    if SHOW_FIGS: plt.show()
    plt.close()

def plot_thermal_with_cumulative(f, Pcorr, P_total_W, title_suffix=""):
    """PSD + cumulative integral in W, also show equivalent <V^2>=R*<P> in title."""
    V2_total = R_OHM * P_total_W
    if f is not None:
        df = np.gradient(f)
        cum_W = np.cumsum(np.clip(Pcorr, 0.0, None) * df)
        x = f / 1e6; xlabel = "Frequency (MHz)"
    else:
        cum_W = np.cumsum(np.clip(Pcorr, 0.0, None))
        x = np.arange(len(Pcorr)); xlabel = "Point index"

    fig, ax = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    ax[0].plot(x, Pcorr, lw=1.6); ax[0].set_ylabel("PSD [W/Hz]")
    ax[0].set_title("Corrected thermal PSD")
    ax[1].plot(x, cum_W, lw=1.6)
    ax[1].set_xlabel(xlabel); ax[1].set_ylabel("Cumulative ∫PSD df [W]")
    ax[1].set_title(f"Total = {P_total_W:.3e} W   (equiv. <V^2> = {V2_total:.3e} V² @ {R_OHM:.0f} Ω)")
    out = os.path.join(OUT_DIR, f"thermal_corrected_cumulative{title_suffix.replace(' ', '_')}.png")
    plt.tight_layout(); plt.savefig(out, dpi=200); print("[saved]", out)
    if SHOW_FIGS: plt.show()
    plt.close()
    return cum_W

def plot_driven_voltage(run_id, y_label="Voltage [V]"):
    f, y, dep, _ = load_trace_1d(run_id)
    plt.figure(figsize=(10, 6))
    if f is None:
        x = np.arange(len(y)); xlabel = "Point index"
    else:
        x = f / 1e6; xlabel = "Frequency (MHz)"
    plt.plot(x, y, "-", lw=1.8, label=f"{dep}")
    plt.xlabel(xlabel); plt.ylabel(y_label); plt.title(f"Driven sweep (run {run_id})")
    plt.legend()
    out_png = os.path.join(OUT_DIR, f"driven_voltage_run{run_id}.png")
    plt.tight_layout(); plt.savefig(out_png, dpi=200); print("[saved]", out_png)
    if SHOW_FIGS: plt.show()
    plt.close()

    # CSV
    if f is None:
        idx = np.arange(len(y), dtype=float)
        save_csv(os.path.join(OUT_DIR, f"driven_run{run_id}.csv"), "index,voltage_V", idx, y)
    else:
        save_csv(os.path.join(OUT_DIR, f"driven_run{run_id}.csv"), "f_Hz,voltage_V", f, y)
    return f, y

# ------------------------ Main ------------------------

def main():
    initialise_or_create_database_at(DB_PATH)

    # Subtraction (kept in W/Hz)
    f, Pdat, Pbg, Pcorr = subtract_background(
        THERMAL_RUN_ID, BACKGROUND_RUN_ID,
        mode=SUBTRACTION_MODE,
        voltage_domain=VOLTAGE_DOMAIN_SUB
    )

    print(f"[thermal] data N={len(Pdat)}, bg N={len(Pbg)}, corrected N={len(Pcorr)}")
    if f is not None and len(f) > 1:
        print(f"[thermal] freq span: {f[0]:.6g} → {f[-1]:.6g} Hz")

    # Plots
    tag = f"(run {THERMAL_RUN_ID} − {BACKGROUND_RUN_ID})"
    plot_thermal_corrected(f, Pcorr, title_suffix=tag)
    plot_thermal_diagnostic(f, Pdat, Pbg, Pcorr, THERMAL_RUN_ID, BACKGROUND_RUN_ID)

    # Integral: <P> in W, then report <V^2> = R * <P>
    if f is None:
        raise RuntimeError("Frequency axis missing: cannot integrate to physical units. Rerun with frequency available.")
    P_total_W = integrate_psd_W_per_Hz(Pcorr, f, one_sided=True, clip_neg=True)
    V2_total  = R_OHM * P_total_W
    print(f"<P>   = {P_total_W:.6e} W")
    print(f"<V^2> = {V2_total:.6e} V^2   (using R = {R_OHM:.0f} Ω)")

    # Cumulative + CSVs
    cum_W = plot_thermal_with_cumulative(f, Pcorr, P_total_W, title_suffix=tag)

    # Save CSVs (all W/Hz)
    save_csv(os.path.join(OUT_DIR, f"thermal_data_run{THERMAL_RUN_ID}.csv"), 
             "f_Hz,Pdata_W_per_Hz", f, Pdat)
    save_csv(os.path.join(OUT_DIR, f"thermal_bg_run{BACKGROUND_RUN_ID}.csv"), 
             "f_Hz,Pbg_W_per_Hz", f, Pbg)
    save_csv(os.path.join(OUT_DIR, f"thermal_corrected_run{THERMAL_RUN_ID}_minus_{BACKGROUND_RUN_ID}.csv"), 
             "f_Hz,Pcorr_W_per_Hz", f, Pcorr)
    save_csv(os.path.join(OUT_DIR, f"thermal_corrected_cumulative_run{THERMAL_RUN_ID}.csv"),
             "f_Hz,cumulative_W", f, cum_W)

    # Driven sweep (voltage)
    if DRIVEN_RUN_ID is not None:
        plot_driven_voltage(DRIVEN_RUN_ID, y_label="Voltage [V]")

    print("Done. All figures and CSVs are in:", OUT_DIR)

if __name__ == "__main__":
    main()

