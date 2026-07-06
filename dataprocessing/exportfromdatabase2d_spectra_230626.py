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

# ============================================================
# Effective mass of mode 1 / mode 2
#   sin(pi x / L) shape  =>  M = M_NT / 2
#   M_NT = (2 * M_C) * (pi * d * L / A_hex)
# ============================================================
M_C   = 12.011 * 1.66053907e-27   # kg, mass of a carbon atom
A_hex = 5.2e-20                   # m^2, area of a graphene hexagon
d_nt  = 2.1e-9                    # m, nanotube diameter
L_nt  = 860e-9                    # m, nanotube length (this device)

M_eff = 0.5 * (2 * M_C) * (np.pi * d_nt * L_nt / A_hex)

# ============================================================
# Equipartition calibration: W/Hz  ->  m^2/Hz
#   <x^2> = (n + 1/2) * hbar / (M * omega_0)
#   integral of measured Lorentzian (W) = A_m * pi * gamma_m
#   calib = <x^2> / integral   [ (m^2/Hz) per (W/Hz) ]
# ============================================================
hbar    = 1.054571817e-34
f_abs   = 121e6                   # Hz, absolute mode frequency
omega_0 = 2 * np.pi * f_abs

n_th    = 8.0                                       # assumed phonon occupancy
x2_mean = (n_th + 0.5) * hbar / (M_eff * omega_0)   # <x^2> in m^2

A_lor_W = A_m * np.pi * gamma_m                     # Lorentzian area [W]
calib   = x2_mean / A_lor_W                         # m^2/Hz per W/Hz

flat_psd_xx  = flat_psd * calib
noise_floor_xx = noise_floor * calib
A_xx   = A_m   * calib
off_xx = off_m * calib

Q = f_abs / fwhm

fit_peak = A_m + off_m              # peak value of the Lorentzian fit
signal = fit_peak - noise_floor     # height above the true noise floor
snr = signal / noise_floor

# ============================================================
# Final plot — displacement PSD in m^2 / Hz
# ============================================================
plt.plot(freq, flat_psd_xx, label="data")
#plt.plot(freq, lorentzian(freq, A_xx, gamma_m, f0_m, off_xx), '--',
#         label=fr"Lorentzian, $\gamma$={gamma_m:.1f} Hz, Q={Q:.0f}")
#plt.axhline(noise_floor_xx, color='k', lw=0.5, ls=':')
plt.xlabel("frequency delta from demod [Hz]")
plt.ylabel(r"$S_{xx}$ [m$^2$/Hz]")
plt.ylim(bottom=0)
#plt.legend(loc="upper right")
plt.show()

print(f"Effective mass M:        {M_eff:.3e} kg")
print(f"Absolute mode freq f0:   {f_abs:.3e} Hz")
print(f"Peak offset from demod:  {f0_m:.2f} Hz")
print(f"Linewidth gamma (HWHM):  {gamma_m:.2f} Hz")
print(f"Q factor:                {Q:.0f}")
print(f"<x^2> (n = {n_th:g}):         {x2_mean:.3e} m^2   (x_rms = {np.sqrt(x2_mean):.3e} m)")
print(f"Lorentzian area:         {A_lor_W:.3e} W")
print(f"Calibration:             {calib:.3e} m^2/W")
print(f"Noise floor:             {noise_floor:.3e} W/Hz  ->  {noise_floor_xx:.3e} m^2/Hz")
print(f"Signal above floor:      {signal:.3e} W/Hz")
print(f"SNR:                     {(0.96*snr):.2f}")