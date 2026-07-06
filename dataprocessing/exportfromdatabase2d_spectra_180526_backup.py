import math
#import math
#import csv
import matplotlib.pyplot as plt
#import os


import pandas as pd
import qcodes as qc
import numpy as np

def voltage_to_dbm(v_rms, rbw, impedance=50):
    # Convert voltage to dBm using the formula
    dbm = 20 * np.log10(v_rms) + 10 * np.log10(1 / (impedance * rbw)) + 30
    return dbm
def voltage_to_psd(v_rms, rbw, impedance=50):

    # Calculate PSD using the formula
    psd = (v_rms ** 2) / (impedance * rbw)
    return psd

def voltage_to_asd(v_rms, rbw):
    # Voltage amplitude spectral density [V/sqrt(Hz)]
    return v_rms / np.sqrt(rbw)

# ---------------- database / run selection ----------------
# Toggle whichever block you want. The older CD11_D7_C1.db runs (≲2615) store
# the FFT data with different setpoint names than the zurichdata.db runs, so
# the extraction below auto-discovers them instead of hard-coding.

# qc.config["core"]["db_location"]=r"D:\databases CD11_D7_C1\CD11_D7_C1_zurichdata.db"
# RUN_ID = 97        # newer-layout dataset

qc.config["core"]["db_location"]=r"D:\databases CD11_D7_C1\CD11_D7_C1.db"
RUN_ID = 2156
#qc.config["core"]["db_location"]=r"D:\databases CD11_D7_C1\CD11_D7_C1_part2.db"
#RUN_ID = 174
       # older-layout dataset (same logic also handles 2614, 2615, ...)

dataset = qc.load_by_id(RUN_ID)
experiments = qc.experiments()
print(f"Loaded run {RUN_ID} from {qc.config['core']['db_location']}")

# ---------------- robust FFT-parameter extraction ----------------
# Find the FFT/voltage parameter and its setpoints without assuming names.
interdeps = dataset.description.interdeps
non_deps = list(interdeps.non_dependencies)
print("Non-dependency parameters:")
for ps in non_deps:
    print(f"  - {ps.name}")

def _looks_like_fft(name):
    n = name.lower()
    return ("fft" in n) or ("voltage" in n) or ("spec" in n)

param_spec = next((ps for ps in non_deps if _looks_like_fft(ps.name)), non_deps[0])
param_name = param_spec.name
print(f"Using parameter: {param_name}")

data_xy = dataset.get_parameter_data(param_spec)
sub = data_xy[param_name]
value_arr = np.asarray(sub[param_name]).ravel()

setpoint_keys = [k for k in sub.keys() if k != param_name]
for k in setpoint_keys:
    a = np.asarray(sub[k]).ravel()
    print(f"  setpoint {k!r}: shape={a.shape}, n_unique={len(np.unique(a))}, "
          f"range=[{a.min():.4g}, {a.max():.4g}]")

def _is_freq(name, arr):
    n = name.lower()
    if any(tok in n for tok in ("freq", "fft")):
        return True
    if any(tok in n for tok in ("time", "rep")):
        return False
    # Fallback heuristic: frequency axes usually have many points spanning many Hz.
    return len(np.unique(arr)) > 10 and float(np.ptp(arr)) > 10.0

freq_key = None
time_key = None
if len(setpoint_keys) == 1:
    freq_key = setpoint_keys[0]
elif len(setpoint_keys) >= 2:
    scored = [(k, _is_freq(k, np.asarray(sub[k]).ravel())) for k in setpoint_keys]
    fc = [k for k, f in scored if f]
    tc = [k for k, f in scored if not f]
    if fc and tc:
        freq_key, time_key = fc[0], tc[0]
    else:
        # Disambiguate by cardinality: the freq axis usually has more unique points.
        ordered = sorted(setpoint_keys,
                         key=lambda k: -len(np.unique(np.asarray(sub[k]).ravel())))
        freq_key = ordered[0]
        time_key = ordered[1] if len(ordered) > 1 else None
else:
    raise RuntimeError(f"No setpoint columns found for {param_name}; cannot reshape.")

print(f"Frequency axis: {freq_key!r}; time axis: {time_key!r}")

freq_arr = np.asarray(sub[freq_key]).ravel()
time_arr = np.asarray(sub[time_key]).ravel() if time_key else None

freq = np.unique(freq_arr)
time = np.unique(time_arr) if time_arr is not None else np.array([0.0])

if time_arr is not None and len(time) > 1:
    nr_time_points = len(time)
else:
    nr_time_points = max(1, round(len(value_arr) / len(freq)))

real_time_axis = np.linspace(0, 4.772 * nr_time_points, nr_time_points)
freq_spectrum_real = np.zeros((nr_time_points, len(freq)))

# Prefer the natural (time-outer, freq-inner) C-order reshape, but only if the
# row layout actually matches. Otherwise scatter values into the grid using the
# setpoint columns directly — this is what the legacy datasets need.
expected_total = nr_time_points * len(freq)
natural_ok = False
if value_arr.size == expected_total:
    try:
        freq_rows = freq_arr.reshape(nr_time_points, len(freq))
        natural_ok = np.allclose(freq_rows, freq_rows[0])
    except Exception:
        natural_ok = False

if natural_ok:
    freq_spectrum_real = value_arr.reshape(nr_time_points, len(freq))
    print("Reshape: natural C-order (time outer, freq inner).")
else:
    print(f"Reshape: scatter-by-index (value_arr={value_arr.size}, grid={expected_total}).")
    freq_to_idx = {float(v): i for i, v in enumerate(freq)}
    if time_arr is not None:
        time_to_idx = {float(v): i for i, v in enumerate(time)}
    for i, v in enumerate(value_arr):
        fi = freq_to_idx[float(freq_arr[i])]
        ti = time_to_idx[float(time_arr[i])] if time_arr is not None else 0
        freq_spectrum_real[ti, fi] = v


plt.pcolor(freq,real_time_axis,freq_spectrum_real)
plt.xlabel("frequency delta from demod [Hz]")
plt.ylabel("time [s]")
plt.show()

average=True
if average:
    time_avg=np.mean(freq_spectrum_real,axis=0)
    plt.plot(freq,time_avg)
    plt.xlabel("frequency delta from demod [Hz]")
    plt.ylabel("Voltage [V]")
    plt.show()

    plt.plot(freq,voltage_to_psd(time_avg,0.209))
    plt.xlabel("frequency delta from demod [Hz]")
    plt.ylabel("Power [W/Hz]")
    #plt.ylim(top=5.5e-13)
    plt.show()

    asd = voltage_to_asd(time_avg, 0.209)
    asd_smoothed = pd.Series(asd).rolling(window=10, center=True, min_periods=1).mean().to_numpy()
    plt.plot(freq, asd_smoothed)
    plt.xlabel("frequency delta from demod [Hz]")
    plt.ylabel(r"Voltage spectral density [V/$\sqrt{\mathrm{Hz}}$] (rolling avg, n=10)")
    plt.show()


from scipy.optimize import curve_fit

def lorentzian(f, A, gamma, f0, offset):
    return A * (gamma**2) / ((f - f0)**2 + gamma**2) + offset


def filter_response(f, A, fc, f0):
    return A / (1 + ((f - f0)/fc)**2)**3


# ---------------- PSD ----------------
psd = voltage_to_psd(time_avg, 0.209)

# ---------------- Background fit (Option A: iterative, 3rd-order cascaded RC) ----------------
peak_halfwidth = 1000   # Hz — initial mask around the mechanical peak
bg_mask = np.abs(freq - freq[np.argmax(psd)]) > peak_halfwidth

for _ in range(5):
    popt_bg, _ = curve_fit(
        filter_response, freq[bg_mask], psd[bg_mask],
        p0=[psd[bg_mask].max(), 2800.0, 0.0], maxfev=10000)
    fit = filter_response(freq, *popt_bg)
    resid = psd - fit
    sigma = np.std(resid[bg_mask])
    new_mask = resid < 2 * sigma     # reject points >2σ above fit
    if np.array_equal(new_mask, bg_mask):
        break
    bg_mask = new_mask

filter_bg = fit
print(f"Background fit: fc = {popt_bg[1]:.1f} Hz, kept {bg_mask.sum()}/{len(freq)} points")

# Diagnostic: data + fitted background
plt.plot(freq, psd, label="data")
plt.plot(freq, filter_bg, '--', label=f"filter fit (n=3, fc={popt_bg[1]:.0f} Hz)")
plt.xlabel("frequency delta from demod [Hz]")
plt.ylabel("Power [W/Hz]")
plt.legend()
plt.show()

# ---------------- Background-subtracted spectrum ----------------
flat_psd = psd - filter_bg

plt.plot(freq, flat_psd)
plt.axhline(0, color='k', lw=0.5)
plt.xlabel("frequency delta from demod [Hz]")
plt.ylabel("Power [W/Hz] (background subtracted)")
plt.show()

# ---------------- Mechanical Lorentzian fit on the flattened spectrum ----------------
# ---------------- Background-subtracted spectrum ----------------
noise_floor = filter_response(0.0, *popt_bg)   # filter value at its center
flat_psd = psd - filter_bg + noise_floor

plt.plot(freq, flat_psd)
plt.axhline(noise_floor, color='k', lw=0.5, ls=':')
plt.xlabel("frequency delta from demod [Hz]")
plt.ylabel("Power [W/Hz] (background subtracted)")
plt.show()

# ---------------- Mechanical Lorentzian fit on the flattened spectrum ----------------
peak_region = np.abs(freq - freq[np.argmax(flat_psd)]) < peak_halfwidth
p0_mech = [flat_psd.max() - noise_floor, 100.0, freq[np.argmax(flat_psd)], noise_floor]
popt_mech, _ = curve_fit(lorentzian, freq[peak_region], flat_psd[peak_region], p0=p0_mech)

A_m, gamma_m, f0_m, off_m = popt_mech
fwhm = 2 * gamma_m

# Set this to your demod LO frequency in Hz to get an absolute Q
f_demod = None
f_abs = f_demod + f0_m if f_demod else None
Q_str = f"Q ≈ {f_abs/fwhm:.0f}" if f_abs else f"Q (needs f_demod) — FWHM = {fwhm:.1f} Hz"

fit_peak = A_m + off_m              # peak value of the Lorentzian fit
signal = fit_peak - noise_floor     # height above the true noise floor
snr = signal / noise_floor

plt.plot(freq, flat_psd, label="background")
plt.plot(freq, 0.96*lorentzian(freq, *popt_mech), '--',
         label=f"Lorentzian fit, γ={gamma_m:.1f} Hz")
plt.xlabel("frequency delta from demod [Hz]")
plt.ylabel("Power [W/Hz]")
plt.ylim(bottom=0)
plt.legend(loc="upper right")
plt.show()

print(f"Peak offset from demod:  {f0_m:.2f} Hz")
print(f"Linewidth γ (HWHM):      {gamma_m:.2f} Hz")
print(f"Lorentzian peak value:   {fit_peak:.3e} W/Hz")
print(f"Noise floor:             {noise_floor:.3e} W/Hz")
print(f"Signal above floor:      {signal:.3e} W/Hz")
print(f"SNR:                     {(0.96*snr):.2f}")