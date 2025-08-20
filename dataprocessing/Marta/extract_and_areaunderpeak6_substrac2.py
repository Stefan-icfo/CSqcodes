import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.interpolate import interp1d

# -----------------------------
# CONFIGURATION: Update these run numbers as needed
# -----------------------------
run_data = 109 # e.g., measurement with drive
run_background = 111 # e.g., background without drive

# -----------------------------
# Database configuration
# -----------------------------
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4.db"
)

# -----------------------------
# Lorentzian function
# -----------------------------
def lorentzian(x, x0, gamma, A, y0):
    return y0 + (A * gamma**2) / ((x - x0)**2 + gamma**2)

# -----------------------------
# Function to extract sorted data
# -----------------------------
def load_sorted_data(run_id, signal_key="avg_avg_psd_nodrive", freq_key="freq_param"):
    dataset = qc.load_by_id(run_id)
    param_data = dataset.get_parameter_data()
    y = param_data[signal_key][signal_key].flatten()
    x = param_data[signal_key][freq_key].flatten() / 1e6  # Hz → MHz
    sorted_indices = np.argsort(x)
    return x[sorted_indices], y[sorted_indices]

# -----------------------------
# Load data
# -----------------------------
x, y = load_sorted_data(run_data)
x_B, y_B = load_sorted_data(run_background)

assert np.allclose(x, x_B), "X-axes do not match between the two runs!"
y_diff = y - y_B

# -----------------------------
# Lorentzian fit
# -----------------------------
x_fit = x
y_fit = y_diff

x0_guess = x_fit[np.argmax(y_fit)]
gamma_guess = 0.00003
A_guess = max(y_fit)
y0_guess = min(y_fit)
p0 = [x0_guess, gamma_guess, A_guess, y0_guess]

popt, pcov = curve_fit(lorentzian, x_fit, y_fit, p0=p0)
perr = np.sqrt(np.diag(pcov))
x0_fit, gamma_fit, A_fit, y0_fit = popt
x0_err, gamma_err, A_err, y0_err = perr

print("\n✅ Fit on subtracted data (±1σ):")
print(f"  x0     = {x0_fit:.6f} ± {x0_err:.2e} MHz")
print(f"  gamma  = {gamma_fit:.6f} ± {gamma_err:.2e} MHz")
print(f"  A      = {A_fit:.3e} ± {A_err:.1e} W/Hz")
print(f"  y0     = {y0_fit:.3e} ± {y0_err:.1e} W/Hz")

# -----------------------------
# Area under fit and offset
# -----------------------------
lorentz_fit = lambda x: lorentzian(x, *popt)
background = lambda x: y0_fit
x1, x2 = x_fit.min(), x_fit.max()

area_total_fit, _ = quad(lorentz_fit, x1, x2)
area_offset_fit, _ = quad(background, x1, x2)
net_area_fit = area_total_fit - area_offset_fit

# -----------------------------
# Area under experimental data
# -----------------------------
y_corr = y_diff - y0_fit
area_data_net = np.trapz(y_corr, x)

print(f"\n📐 Area under Lorentzian fit:")
print(f"  Total area (with offset):     {area_total_fit:.3e} W/Hz·MHz")
print(f"  Offset area (flat y0):        {area_offset_fit:.3e} W/Hz·MHz")
print(f"  ✅ Net peak area (fit):        {net_area_fit:.3e} W/Hz·MHz")

print(f"\n📐 Area under subtracted data:")
print(f"  ✅ Net peak area (data - y0):  {area_data_net:.3e} W/Hz·MHz")

# -----------------------------
# Monte Carlo estimate of fit area
# -----------------------------
N_samples = 1000
rng = np.random.default_rng(seed=42)
samples = rng.multivariate_normal(popt, pcov, size=N_samples)
area_samples = [
    quad(lambda x: lorentzian(x, *params), x1, x2)[0] - quad(lambda x: params[3], x1, x2)[0]
    for params in samples
]
area_samples = np.array(area_samples)
area_mean = np.mean(area_samples)
area_std = np.std(area_samples)

print(f"\n📊 Monte Carlo estimate of fit area:")
print(f"  Mean net area:     {area_mean:.3e} W/Hz·MHz")
print(f"  Std deviation:     {area_std:.3e} W/Hz·MHz")

# -----------------------------
# PLOT 1: Fit + area under Lorentzian
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x, y_diff, 'o', color='green', markersize=3, label=f'Subtracted ({run_data} - {run_background})')

x_dense = np.linspace(x1, x2, 1000)
y_lorentz = lorentz_fit(x_dense)
y_offset = np.full_like(x_dense, y0_fit)

plt.plot(x_dense, y_lorentz, '-', color='red', linewidth=2, label="Lorentzian fit")
plt.fill_between(x_dense, y_lorentz, color='lightgreen', alpha=0.4, label='Total area (fit)')
plt.fill_between(x_dense, y_offset, y_lorentz, where=(y_lorentz > y_offset),
                 color='yellow', alpha=0.6, label='Net peak area (fit - offset)')

fit_text = (
    f"x₀ = {x0_fit:.4f} ± {x0_err:.1e} MHz\n"
    f"γ = {gamma_fit:.5f} ± {gamma_err:.1e} MHz\n"
    f"A = {A_fit:.2e} ± {A_err:.1e} W/Hz\n"
    f"y₀ = {y0_fit:.2e} ± {y0_err:.1e} W/Hz\n"
    f"Area (fit) = {net_area_fit:.2e}\n"
    f"Area (data) = {area_data_net:.2e} W/Hz·MHz"
)
plt.text(0.02, 0.95, fit_text, transform=plt.gca().transAxes,
         fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.xlabel("Frequency (MHz)", fontsize=14)
plt.ylabel("PSD (W/Hz)", fontsize=14)
plt.title(f"Lorentzian fit on subtracted spectrum (Run {run_data})", fontsize=16)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# -----------------------------
# PLOT 2: Experimental subtracted data with area (blue fill)
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x, y_corr, color='blue', linewidth=2, label='Subtracted data - offset')
plt.fill_between(x, 0, y_corr, where=(y_corr > 0), color='skyblue', alpha=0.6, label='Net area (data)')

plt.xlabel("Frequency (MHz)", fontsize=14)
plt.ylabel("PSD - offset (W/Hz)", fontsize=14)
plt.title(f"Area under subtracted data (Run {run_data} - {run_background})", fontsize=16)
plt.text(0.02, 0.95, f"Area = {area_data_net:.2e} W/Hz·MHz", transform=plt.gca().transAxes,
         fontsize=11, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
