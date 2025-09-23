#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.interpolate import interp1d

# -----------------------------
# USER CONFIG
# -----------------------------
run_data = 1480        # measurement with drive
run_background = 1465  # background without drive
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v8.db"

# QCoDeS keys used by your dataset
signal_key = "avg_avg_psd_nodrive"
freq_key   = "freq_param"

# Optional: restrict fit to a window around the peak (MHz); set None to fit all
fit_window_MHz = None  # e.g. (94.7, 95.1)

# -----------------------------
# Setup DB
# -----------------------------
qc.config["core"]["db_location"] = db_path

# -----------------------------
# Model
# -----------------------------
def lorentzian(x, x0, gamma, A, y0):
    """y0 + A*gamma^2 / ((x-x0)^2 + gamma^2)"""
    return y0 + (A * gamma**2) / ((x - x0)**2 + gamma**2)

# -----------------------------
# IO helpers
# -----------------------------
def load_sorted_data(run_id, signal_key=signal_key, freq_key=freq_key):
    ds = qc.load_by_id(run_id)
    pd = ds.get_parameter_data()
    y = np.asarray(pd[signal_key][signal_key]).ravel()
    x = np.asarray(pd[signal_key][freq_key]).ravel() / 1e6  # Hz -> MHz
    idx = np.argsort(x)
    return x[idx], y[idx]

# -----------------------------
# Load
# -----------------------------
x,   y   = load_sorted_data(run_data)
x_B, y_B = load_sorted_data(run_background)

# -----------------------------
# Interpolate background onto data grid
# (handles different frequency grids gracefully)
# -----------------------------
# Use linear interpolation; outside the background range we keep NaN so those
# points are ignored in subtraction/fit.
interpB = interp1d(
    x_B, y_B, kind="linear",
    bounds_error=False, fill_value=np.nan, assume_sorted=True
)
y_B_on_x = interpB(x)

# Make subtraction only where both are finite
mask = np.isfinite(y) & np.isfinite(y_B_on_x)
x = x[mask]
y = y[mask]
y_diff = y - y_B_on_x[mask]

# -----------------------------
# Optional: restrict to fit window
# -----------------------------
if fit_window_MHz is not None:
    a, b = fit_window_MHz
    mfit = (x >= a) & (x <= b)
    x_fit = x[mfit]
    y_fit = y_diff[mfit]
else:
    x_fit = x
    y_fit = y_diff

if x_fit.size < 5:
    raise RuntimeError("Not enough points in the (optional) fit window to perform a fit.")

# -----------------------------
# Lorentzian fit (with bounds)
# -----------------------------
x0_guess = x_fit[np.argmax(y_fit)]
gamma_guess = max(1e-5, 0.00003)   # MHz; positive
A_guess = np.nanmax(y_fit) if np.isfinite(np.nanmax(y_fit)) else 1.0
y0_guess = np.nanmedian(y_fit)

p0 = [x0_guess, gamma_guess, A_guess, y0_guess]
bounds = (
    [x_fit.min(),  0.0,       -np.inf, -np.inf],   # lower: gamma >= 0
    [x_fit.max(),  np.inf,     np.inf,   np.inf],   # upper
)

popt, pcov = curve_fit(lorentzian, x_fit, y_fit, p0=p0, bounds=bounds, maxfev=20000)
perr = np.sqrt(np.diag(pcov))
x0_fit, gamma_fit, A_fit, y0_fit = popt
x0_err, gamma_err, A_err, y0_err = perr

print("\n✅ Fit on subtracted data (±1σ):")
print(f"  x0     = {x0_fit:.6f} ± {x0_err:.2e} MHz")
print(f"  gamma  = {gamma_fit:.6f} ± {gamma_err:.2e} MHz")
print(f"  A      = {A_fit:.3e} ± {A_err:.1e} W/Hz")
print(f"  y0     = {y0_fit:.3e} ± {y0_err:.1e} W/Hz")

# -----------------------------
# Areas
# -----------------------------
lorentz_fit = lambda xx: lorentzian(xx, *popt)
x1, x2 = (x_fit.min(), x_fit.max())

# Fit area
area_total_fit, _  = quad(lorentz_fit, x1, x2)
area_offset_fit, _ = quad(lambda xx: y0_fit, x1, x2)
net_area_fit = area_total_fit - area_offset_fit

# Data area (subtract fitted offset y0 only)
y_corr = y_diff - y0_fit
area_data_net = np.trapz(y_corr[(x >= x1) & (x <= x2)], x[(x >= x1) & (x <= x2)])

print(f"\n📐 Areas over [{x1:.6f}, {x2:.6f}] MHz:")
print(f"  Net peak area (fit):  {net_area_fit:.3e} W/Hz·MHz")
print(f"  Net peak area (data): {area_data_net:.3e} W/Hz·MHz")

# -----------------------------
# Monte Carlo uncertainty of fit area
# -----------------------------
N_samples = 1000
rng = np.random.default_rng(42)
samples = rng.multivariate_normal(popt, pcov, size=N_samples)
areas = []
for p in samples:
    f = lambda xx: lorentzian(xx, *p)
    at, _ = quad(f, x1, x2)
    ao, _ = quad(lambda xx: p[3], x1, x2)  # y0 from sample
    areas.append(at - ao)
areas = np.array(areas)
print(f"\n📊 Monte Carlo net area: {areas.mean():.3e} ± {areas.std():.3e} W/Hz·MHz (1σ)")

# -----------------------------
# PLOT 1: Subtraction & fit
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x, y_diff, 'o', ms=3, alpha=0.8, label=f"Subtracted ({run_data} − {run_background})")

xd = np.linspace(x1, x2, 1500)
yl = lorentz_fit(xd)
yoff = np.full_like(xd, y0_fit)

plt.plot(xd, yl, '-', lw=2, color='crimson', label='Lorentzian fit')
plt.fill_between(xd, yl, color='lightgreen', alpha=0.35, label='Total area (fit)')
plt.fill_between(xd, yoff, yl, where=(yl > yoff), color='gold', alpha=0.55, label='Net area (fit − offset)')

txt = (
    f"x₀ = {x0_fit:.5f} ± {x0_err:.1e} MHz\n"
    f"γ  = {gamma_fit:.5f} ± {gamma_err:.1e} MHz\n"
    f"A  = {A_fit:.2e} ± {A_err:.1e} W/Hz\n"
    f"y₀ = {y0_fit:.2e} ± {y0_err:.1e} W/Hz\n"
    f"Net area (fit)  = {net_area_fit:.2e}\n"
    f"Net area (data) = {area_data_net:.2e} W/Hz·MHz"
)
plt.text(0.02, 0.98, txt, transform=plt.gca().transAxes, va='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.85), fontsize=10)

plt.xlabel("Frequency (MHz)")
plt.ylabel("PSD (W/Hz)")
plt.title(f"Lorentzian fit on subtracted spectrum (run {run_data} − {run_background})")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# -----------------------------
# PLOT 2: Net data area (data − fitted offset)
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x, y_corr, lw=2, color='navy', label='Subtracted data − fitted offset')
plt.fill_between(x, 0, y_corr, where=(y_corr > 0), color='skyblue', alpha=0.6, label='Net area (data)')
plt.xlabel("Frequency (MHz)")
plt.ylabel("PSD − y₀ (W/Hz)")
plt.title(f"Area under subtracted data (run {run_data} − {run_background})")
plt.text(0.02, 0.97, f"Area = {area_data_net:.2e} W/Hz·MHz",
         transform=plt.gca().transAxes, va='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.85), fontsize=11)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
