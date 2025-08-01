import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit
from scipy.integrate import quad
from numpy.random import default_rng

# -----------------------------
# CONFIGURATION
# -----------------------------
run_data = 109
run_background = 113

# -----------------------------
# Database location
# -----------------------------
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v2.db"
)

# -----------------------------
# Lorentzian model
# -----------------------------
def lorentzian(x, x0, gamma, A, y0):
    return y0 + (A * gamma**2) / ((x - x0)**2 + gamma**2)

# -----------------------------
# Load and sort data
# -----------------------------
def load_sorted_data(run_id, signal_key="avg_avg_psd_nodrive", freq_key="freq_param"):
    ds = qc.load_by_id(run_id)
    pdata = ds.get_parameter_data()
    y = pdata[signal_key][signal_key].flatten()
    x = pdata[signal_key][freq_key].flatten() / 1e6  # Convert Hz to MHz
    idx = np.argsort(x)
    return x[idx], y[idx]

x, P_total = load_sorted_data(run_data)
x_bkg, P_noise = load_sorted_data(run_background)
assert np.allclose(x, x_bkg), "X axes mismatch!"

# -----------------------------
# Correct noise subtraction (voltage domain)
# -----------------------------
V_total = np.sqrt(P_total)
V_noise = np.sqrt(P_noise)
V_signal = V_total - V_noise
P_signal = V_signal**2

# -----------------------------
# Fit Lorentzian
# -----------------------------
x_fit = x
y_fit = P_signal
x0_guess = x_fit[np.argmax(y_fit)]
p0 = [x0_guess, 0.00003, max(y_fit), min(y_fit)]
popt, pcov = curve_fit(lorentzian, x_fit, y_fit, p0=p0)
perr = np.sqrt(np.diag(pcov))
x0_fit, gamma_fit, A_fit, y0_fit = popt
x0_err, gamma_err, A_err, y0_err = perr

print("\n✅ Lorentzian Fit (±1σ):")
print(f"  x₀     = {x0_fit:.6f} ± {x0_err:.2e} MHz")
print(f"  γ      = {gamma_fit:.6f} ± {gamma_err:.2e} MHz")
print(f"  A      = {A_fit:.2e} ± {A_err:.1e} W/Hz")
print(f"  y₀     = {y0_fit:.2e} ± {y0_err:.1e} W/Hz")

# -----------------------------
# Compute areas
# -----------------------------
x1, x2 = x_fit.min(), x_fit.max()
lorentz_fit = lambda x: lorentzian(x, *popt)
background = lambda x: y0_fit
area_total, _ = quad(lorentz_fit, x1, x2)
area_offset, _ = quad(background, x1, x2)
net_area_fit = area_total - area_offset
y_corr = y_fit - y0_fit
area_data_net = np.trapz(y_corr, x)

print(f"\n📐 Area under Lorentzian fit:")
print(f"  Total area:       {area_total:.3e} W/Hz·MHz")
print(f"  Offset area:      {area_offset:.3e} W/Hz·MHz")
print(f"  ✅ Net peak area:  {net_area_fit:.3e} W/Hz·MHz")

print(f"\n📐 Area from data (P_signal - y₀):")
print(f"  ✅ Net area:       {area_data_net:.3e} W/Hz·MHz")

# -----------------------------
# Monte Carlo uncertainty (fit area)
# -----------------------------
rng = default_rng(seed=42)
samples = rng.multivariate_normal(popt, pcov, size=1000)
area_samples = [
    quad(lambda x: lorentzian(x, *p), x1, x2)[0] - quad(lambda x: p[3], x1, x2)[0]
    for p in samples
]
area_mean = np.mean(area_samples)
area_std = np.std(area_samples)

print(f"\n📊 Monte Carlo estimate of area:")
print(f"  Mean net area:    {area_mean:.3e} W/Hz·MHz")
print(f"  Std deviation:    {area_std:.3e} W/Hz·MHz")

# -----------------------------
# Plot 1: Fit + Areas
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x, P_signal, 'o', color='darkgreen', label='Corrected spectrum', markersize=3)

x_dense = np.linspace(x1, x2, 1000)
y_fit_dense = lorentz_fit(x_dense)
plt.plot(x_dense, y_fit_dense, '-', color='red', label='Lorentzian fit')
plt.fill_between(x_dense, y_fit_dense, color='lightgreen', alpha=0.3, label='Total area')
plt.fill_between(x_dense, y0_fit, y_fit_dense, where=(y_fit_dense > y0_fit),
                 color='orange', alpha=0.5, label='Net area')

info = (
    f"x₀ = {x0_fit:.4f} ± {x0_err:.1e} MHz\n"
    f"γ = {gamma_fit:.5f} ± {gamma_err:.1e} MHz\n"
    f"A = {A_fit:.2e} ± {A_err:.1e} W/Hz\n"
    f"y₀ = {y0_fit:.2e} ± {y0_err:.1e} W/Hz\n"
    f"Net area = {net_area_fit:.2e} W/Hz·MHz"
)
plt.text(0.02, 0.95, info, transform=plt.gca().transAxes, fontsize=10,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.xlabel("Frequency (MHz)")
plt.ylabel("PSD (W/Hz)")
plt.title(f"Lorentzian Fit after Voltage-Domain Subtraction (Run {run_data})")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# -----------------------------
# Plot 2: Area under corrected data
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x, y_corr, label='P_signal - y₀', color='blue')
plt.fill_between(x, 0, y_corr, where=(y_corr > 0), color='skyblue', alpha=0.5, label='Net area')
plt.text(0.02, 0.95, f"Net area = {area_data_net:.2e} W/Hz·MHz",
         transform=plt.gca().transAxes, fontsize=11, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.xlabel("Frequency (MHz)")
plt.ylabel("Corrected PSD (W/Hz)")
plt.title("Area under corrected PSD (voltage-based subtraction)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show() 