#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step 1: subtract thermal background and plot corrected PSD
Step 2: integrate the corrected PSD to get <V^2> (or <P>) and show a cumulative integral
Step 3: plot driven sweep (voltage vs frequency)

NO unit conversion is performed: spectra are used exactly as stored in the DB.
You choose the axis/legend units via PSD_UNITS below so figures/CSV headers are correct.
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

# Background subtraction mode:
SUBTRACTION_MODE   = "interpolate"  # "index" (point-by-point) or "interpolate" (recommended)
VOLTAGE_DOMAIN_SUB = False          # False: P_data - P_bg (keeps units as-is)
                                    # True:  (sqrt(P_data) - sqrt(P_bg))^2  (assumes PSD ∝ V^2)

# Plot x-axis: use point index (if frequency missing) or frequency when available
USE_INDEX_AXIS_FOR_PLOT = False

# Show figures on screen?
SHOW_FIGS = True

# Output folder
OUT_DIR = os.path.abspath("./fig_subtract_and_driven_with_integral")
os.makedirs(OUT_DIR, exist_ok=True)
print("Output folder:", OUT_DIR)
try:
    os.startfile(OUT_DIR)  # Windows convenience; ignored elsewhere
except Exception:
    pass

# (Optional) Force the dependent key to read (skip auto-detection), e.g.:
# FORCE_DEP_KEY = "avg_avg_psd_nodrive"
FORCE_DEP_KEY = None
# (Optional) Force the frequency setpoint key (skip auto-detection)
FORCE_FREQ_KEY = None

# ======= Units to display in plots/CSVs (no conversion is applied) =======
# Set this to match what is stored in your DB:
#   "V^2/Hz"  -> integral is "<V^2> [V^2]"
#   "W/Hz"    -> integral is "<P> [W]"
PSD_UNITS = "V^2/Hz"   # or "W/Hz"
# ========================================================================

mpl.rcParams["font.family"] = "Calibri"

# Keys commonly found in these datasets (adjust if needed)
CANDIDATE_SIGNAL_KEYS = [
    "avg_avg_psd_nodrive", "avg_psd_nodrive",
    "avg_avg_psd", "avg_psd",
    "V_fft_avg_avg", "PSD", "psd", "Svv", "V_r", "Amplitude"
]
FREQ_KEY_CANDIDATES = ["freq_param", "compressed_freq_real", "freq", "frequency", "Frequency", "f"]


# ---------------------- Helpers ----------------------

def integral_units_from_psd_units(psd_units: str) -> str:
    """Return the correct units for the integral based on PSD units."""
    u = psd_units.strip().lower().replace(" ", "")
    if u == "v^2/hz" or u == "v2/hz":
        return "V^2"
    if u == "w/hz":
        return "W"
    # Fallback: write as "(PSD units)·Hz"
    return f"{psd_units}·Hz"

def save_csv(path, header, *cols):
    """Save columns as CSV with a header."""
    arr = np.column_stack(cols)
    np.savetxt(path, arr, delimiter=",", header=header, comments="")
    print("[saved]", path)

def load_trace_1d(run_id):
    """
    Load one dependent 1D trace and (if present) its frequency setpoint.
    Keeps acquisition ORDER (no sorting).
    Returns (f_Hz or None, y_values, dep_key, freq_key).
    """
    ds = load_by_id(run_id)
    pd = ds.get_parameter_data()

    # Choose dependent
    if FORCE_DEP_KEY is not None:
        dep = FORCE_DEP_KEY
        if dep not in pd:
            raise KeyError(f"Forced dep key '{dep}' not found in dataset.")
    else:
        dep = next((k for k in CANDIDATE_SIGNAL_KEYS if k in pd), next(iter(pd.keys())))
    blk = pd[dep]

    # Choose frequency setpoint in the same block
    if FORCE_FREQ_KEY is not None:
        freq_key = FORCE_FREQ_KEY
        if freq_key not in blk:
            raise KeyError(f"Forced freq key '{freq_key}' not found in block for dep '{dep}'.")
    else:
        freq_key = None
        for k in blk.keys():
            if k == dep:
                continue
            if any(tag.lower() in k.lower() for tag in FREQ_KEY_CANDIDATES):
                freq_key = k
                break

    y = np.asarray(blk[dep]).ravel().astype(float)
    f = np.asarray(blk[freq_key]).ravel().astype(float) if freq_key is not None else None

    # Drop NaNs while keeping alignment
    m = np.isfinite(y)
    if f is not None:
        m &= np.isfinite(f)
        f = f[m]
    y = y[m]

    print(f"[run {run_id}] dep='{dep}', fkey='{freq_key}', N={len(y)}")
    return f, y, dep, freq_key

def subtract_background(thermal_id, bg_id, mode="interpolate", voltage_domain=False):
    """
    Return aligned axis and PSDs (no conversion):
      f_used (Hz or None), P_data, P_bg_aligned, P_corrected
    - mode='index': truncate to min length and subtract point-by-point (acquisition order)
    - mode='interpolate': interpolate background onto the data frequency grid (recommended)
    - voltage_domain=False: P_corr = P_data - P_bg
      voltage_domain=True : P_corr = (sqrt(P_data) - sqrt(P_bg))^2
    """
    f_d, y_d, _, _ = load_trace_1d(thermal_id)
    f_b, y_b, _, _ = load_trace_1d(bg_id)

    P_d = np.asarray(y_d, dtype=float)
    P_b = np.asarray(y_b, dtype=float)

    if mode.lower().startswith("index"):
        n = min(len(P_d), len(P_b))
        Pd = P_d[:n]
        Pb = P_b[:n]
        if voltage_domain:
            Vd = np.sqrt(np.clip(Pd, 0, None))
            Vb = np.sqrt(np.clip(Pb, 0, None))
            Pc = (Vd - Vb) ** 2
        else:
            Pc = Pd - Pb
        f_used = (None if (f_d is None or USE_INDEX_AXIS_FOR_PLOT) else f_d[:n])
        return f_used, Pd[:n], Pb[:n], Pc

    # Interpolate background onto the data frequency axis
    if f_d is None or f_b is None:
        raise ValueError("Mode 'interpolate' requires both frequency axes.")
    Pbg_interp = np.interp(f_d, f_b, P_b, left=P_b[0], right=P_b[-1])
    if voltage_domain:
        Vd = np.sqrt(np.clip(P_d, 0, None))
        Vb = np.sqrt(np.clip(Pbg_interp, 0, None))
        Pc = (Vd - Vb) ** 2
    else:
        Pc = P_d - Pbg_interp
    return f_d, P_d, Pbg_interp, Pc

def trapz_integral(y, x=None):
    """Robust trapezoidal integration with NaN masking."""
    y = np.asarray(y, dtype=float)
    if x is None:
        # Integrate w.r.t. index (Δ=1) if frequency is missing
        m = np.isfinite(y)
        return np.trapz(y[m])
    x = np.asarray(x, dtype=float)
    m = np.isfinite(y) & np.isfinite(x)
    if np.count_nonzero(m) < 2:
        return np.nan
    return np.trapz(y[m], x[m])

def plot_thermal_corrected(f, P_corr, title_suffix=""):
    """Plot ONLY the background-subtracted thermal PSD."""
    plt.figure(figsize=(10, 6))
    if f is None or USE_INDEX_AXIS_FOR_PLOT:
        x = np.arange(len(P_corr)); xlabel = "Point index"
    else:
        x = f / 1e6; xlabel = "Frequency (MHz)"
    plt.plot(x, P_corr, "-", lw=1.8, label="Thermal (BG-subtracted)")
    plt.xlabel(xlabel)
    plt.ylabel(f"PSD [{PSD_UNITS}]")
    plt.title(f"Thermal PSD after BG subtraction {title_suffix}")
    plt.legend()
    out = os.path.join(OUT_DIR, f"thermal_corrected{title_suffix.replace(' ', '_')}.png")
    plt.tight_layout(); plt.savefig(out, dpi=200)
    print("[saved]", out)
    if SHOW_FIGS:
        plt.show()
    plt.close()

def plot_thermal_diagnostic(f, Pdat, Pbg, Pcorr, thermal_id, bg_id):
    """Diagnostic plot: data vs background vs corrected, to sanity-check subtraction."""
    plt.figure(figsize=(10, 6))
    if f is None or USE_INDEX_AXIS_FOR_PLOT:
        n = len(Pcorr)
        x = np.arange(n); xlabel = "Point index"
        plt.plot(x, Pdat[:n], alpha=0.7, label="Data")
        plt.plot(x, Pbg[:n],  alpha=0.7, label="Background (aligned)")
        plt.plot(x, Pcorr, lw=1.8, label="Corrected")
    else:
        x = f / 1e6; xlabel = "Frequency (MHz)"
        plt.plot(x, Pdat, alpha=0.7, label="Data")
        plt.plot(x, Pbg,  alpha=0.7, label="Background (aligned)")
        plt.plot(x, Pcorr, lw=1.8, label="Corrected")
    plt.xlabel(xlabel)
    plt.ylabel(f"PSD [{PSD_UNITS}]")
    plt.title("Thermal PSD: data vs background vs corrected")
    plt.legend()
    out = os.path.join(OUT_DIR, f"thermal_data_bg_corrected_{thermal_id}_{bg_id}.png")
    plt.tight_layout(); plt.savefig(out, dpi=200)
    print("[saved]", out)
    if SHOW_FIGS:
        plt.show()
    plt.close()

def plot_thermal_with_cumulative(f, Pcorr, title_suffix=""):
    """PSD (top) + cumulative integral (bottom), both with correct units in labels."""
    INT_UNITS = integral_units_from_psd_units(PSD_UNITS)

    # Build cumulative integral vs frequency (or index if f is None)
    if f is not None:
        df = np.gradient(f)
        # Clip negative corrections before integrating (optional, but avoids small negative tails)
        val = np.clip(Pcorr, 0.0, None) * df
        cumulative = np.cumsum(val)
        x_plot = f / 1e6
        xlabel = "Frequency (MHz)"
    else:
        val = np.clip(Pcorr, 0.0, None)  # Δ=1 index
        cumulative = np.cumsum(val)
        x_plot = np.arange(len(Pcorr))
        xlabel = "Point index"

    total_area = float(cumulative[-1]) if len(cumulative) else np.nan

    fig, ax = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    ax[0].plot(x_plot, Pcorr, lw=1.6)
    ax[0].set_ylabel(f"PSD [{PSD_UNITS}]")
    ax[0].set_title("Corrected thermal PSD")

    ax[1].plot(x_plot, cumulative, lw=1.6)
    ax[1].set_xlabel(xlabel)
    ax[1].set_ylabel(f"Cumulative integral [{INT_UNITS}]")
    ax[1].set_title(f"Cumulative ∫ PSD df   (Total = {total_area:.3e} {INT_UNITS})")

    out = os.path.join(OUT_DIR, f"thermal_corrected_cumulative{title_suffix.replace(' ', '_')}.png")
    plt.tight_layout(); plt.savefig(out, dpi=200)
    print("[saved]", out)
    if SHOW_FIGS:
        plt.show()
    plt.close()

    return total_area, cumulative

def plot_driven_voltage(run_id, y_label="Voltage [V]"):
    """
    Plot the driven sweep as stored (no calibration).
    Also saves a CSV and returns (f_Hz or None, y_raw).
    """
    f, y, dep, fkey = load_trace_1d(run_id)

    plt.figure(figsize=(10, 6))
    if f is None:
        x = np.arange(len(y)); xlabel = "Point index"
    else:
        x = f / 1e6; xlabel = "Frequency (MHz)"
    plt.plot(x, y, "-", lw=1.8, label=f"{dep}")
    plt.xlabel(xlabel)
    plt.ylabel(y_label)
    plt.title(f"Driven sweep (run {run_id})")
    plt.legend()
    out_png = os.path.join(OUT_DIR, f"driven_voltage_run{run_id}.png")
    plt.tight_layout(); plt.savefig(out_png, dpi=200)
    print("[saved]", out_png)
    if SHOW_FIGS:
        plt.show()
    plt.close()

    # Save CSV
    if f is None:
        idx = np.arange(len(y), dtype=float)
        save_csv(os.path.join(OUT_DIR, f"driven_run{run_id}.csv"),
                 "index,voltage", idx, y)
    else:
        save_csv(os.path.join(OUT_DIR, f"driven_run{run_id}.csv"),
                 "f_Hz,voltage", f, y)

    # Console ranges
    try:
        if f is not None and len(f) > 0:
            print(f"[driven] f range: {f[0]:.6g} Hz → {f[-1]:.6g} Hz, N={len(f)}")
        print(f"[driven] voltage range: {np.nanmin(y):.3e} → {np.nanmax(y):.3e} V (as stored)")
    except Exception:
        pass

    return f, y


# ------------------------ Main ------------------------

def main():
    initialise_or_create_database_at(DB_PATH)

    # Step 1: background subtraction (no conversion)
    try:
        f, Pdat, Pbg, Pcorr = subtract_background(
            THERMAL_RUN_ID, BACKGROUND_RUN_ID,
            mode=SUBTRACTION_MODE,
            voltage_domain=VOLTAGE_DOMAIN_SUB
        )
    except Exception as e:
        print("Error during background subtraction:", repr(e))
        return

    # Console sanity info
    print(f"[thermal] data N={len(Pdat)}, bg N={len(Pbg)}, corrected N={len(Pcorr)}")
    try:
        print(f"[thermal] corrected range: [{np.nanmin(Pcorr):.3e}, {np.nanmax(Pcorr):.3e}] [{PSD_UNITS}]")
        if f is not None and len(f) > 0:
            print(f"[thermal] frequency span: {f[0]:.6g} Hz → {f[-1]:.6g} Hz (N={len(f)})")
    except Exception:
        pass

    # Step 2: plots and integral
    plot_thermal_corrected(f, Pcorr, title_suffix=f"(run {THERMAL_RUN_ID} − {BACKGROUND_RUN_ID})")
    plot_thermal_diagnostic(f, Pdat, Pbg, Pcorr, THERMAL_RUN_ID, BACKGROUND_RUN_ID)

    # Compute the integral over the corrected PSD (clip negative small tails to 0)
    INT_UNITS = integral_units_from_psd_units(PSD_UNITS)
    if f is not None:
        total_area = trapz_integral(np.clip(Pcorr, 0.0, None), x=f)
    else:
        total_area = trapz_integral(np.clip(Pcorr, 0.0, None), x=None)  # integrates vs index

    print(f"[thermal] Integral of corrected PSD = {total_area:.6e} {INT_UNITS}")

    # Cumulative integral plot + CSV
    total_area2, cumulative = plot_thermal_with_cumulative(f, Pcorr,
                                                           title_suffix=f"(run {THERMAL_RUN_ID} − {BACKGROUND_RUN_ID})")
    # Save CSVs
    try:
        if f is None or USE_INDEX_AXIS_FOR_PLOT:
            idx = np.arange(len(Pcorr), dtype=float)
            save_csv(os.path.join(OUT_DIR, f"thermal_data_run{THERMAL_RUN_ID}.csv"),
                     f"index,Pdata_{PSD_UNITS}", idx, Pdat[:len(Pcorr)])
            save_csv(os.path.join(OUT_DIR, f"thermal_bg_run{BACKGROUND_RUN_ID}.csv"),
                     f"index,Pbg_{PSD_UNITS}", idx, Pbg[:len(Pcorr)])
            save_csv(os.path.join(OUT_DIR, f"thermal_corrected_run{THERMAL_RUN_ID}_minus_{BACKGROUND_RUN_ID}.csv"),
                     f"index,Pcorr_{PSD_UNITS}", idx, Pcorr)
            save_csv(os.path.join(OUT_DIR, f"thermal_corrected_cumulative_run{THERMAL_RUN_ID}.csv"),
                     f"index,cumulative_{INT_UNITS}", idx, cumulative)
        else:
            save_csv(os.path.join(OUT_DIR, f"thermal_data_run{THERMAL_RUN_ID}.csv"),
                     f"f_Hz,Pdata_{PSD_UNITS}", f, Pdat)
            save_csv(os.path.join(OUT_DIR, f"thermal_bg_run{BACKGROUND_RUN_ID}.csv"),
                     f"f_Hz,Pbg_{PSD_UNITS}", f, Pbg)
            save_csv(os.path.join(OUT_DIR, f"thermal_corrected_run{THERMAL_RUN_ID}_minus_{BACKGROUND_RUN_ID}.csv"),
                     f"f_Hz,Pcorr_{PSD_UNITS}", f, Pcorr)
            save_csv(os.path.join(OUT_DIR, f"thermal_corrected_cumulative_run{THERMAL_RUN_ID}.csv"),
                     f"f_Hz,cumulative_{INT_UNITS}", f, cumulative)
    except Exception as e:
        print("CSV save error:", repr(e))

    # Step 3: plot driven sweep (as stored)
    if DRIVEN_RUN_ID is not None:
        try:
            _ = plot_driven_voltage(DRIVEN_RUN_ID, y_label="Voltage [V]")
        except Exception as e:
            print("Error plotting driven sweep:", repr(e))
    else:
        print("DRIVEN_RUN_ID=None -> skipping driven plot.")

    print(f"Done. Integral of corrected PSD = {total_area:.6e} {INT_UNITS}")
    print("All figures and CSVs are in:", OUT_DIR)

if __name__ == "__main__":
    main()
