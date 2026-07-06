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

#database location
qc.config["core"]["db_location"]=r"D:\databases CD11_D7_C1\CD11_D7_C1_zurichdata.db"
experiments=qc.experiments()

dataset=qc.load_by_id(97)


pdf=dataset.to_pandas_dataframe_dict()
freq_spec_raw=pdf["Voltage_fft_avg"]
freq_spec_np=np.array(freq_spec_raw)
# ---------------------Geting the data from the database---------------------
# pprint(dataset.get_parameter_data())
interdeps = dataset.description.interdeps
param_spec = interdeps.non_dependencies[0]  # hall resistance data
param_name = param_spec.name
data_xy = dataset.get_parameter_data(param_spec)
xy = data_xy[param_name][param_name]

#g1:outer gate
#g2:inner gate

time_raw = data_xy[param_name]['time_param']
freq_raw = data_xy[param_name]['freq_param']
time_np=np.array(time_raw)
freq_np=np.array(freq_raw)
time=np.unique(time_np)
freq=np.unique(freq_np)
#real_time_len
nr_time_points=round(len(freq_spec_np)/len(freq))
real_time_axis=np.linspace(0,4.772*nr_time_points,nr_time_points)

freq_spectrum_real=np.zeros([nr_time_points, len(freq)])


for m in range(nr_time_points):
    for n in range(len(freq)):
        freq_spectrum_real[m,n]=freq_spec_np[m*len(freq)+n]

        #integral


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
# plt.plot(freq, 0.96*lorentzian(freq, *popt_mech), '--',
         #label=f"Lorentzian fit, γ={gamma_m:.1f} Hz")
plt.xlabel("frequency delta from demod [Hz]")
plt.ylabel("Power [W/Hz]")
plt.ylim(bottom=0)
#plt.legend(loc="upper right")
plt.show()

print(f"Peak offset from demod:  {f0_m:.2f} Hz")
print(f"Linewidth γ (HWHM):      {gamma_m:.2f} Hz")
print(f"Lorentzian peak value:   {fit_peak:.3e} W/Hz")
print(f"Noise floor:             {noise_floor:.3e} W/Hz")
print(f"Signal above floor:      {signal:.3e} W/Hz")
print(f"SNR:                     {(0.96*snr):.2f}")