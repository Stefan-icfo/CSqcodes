import numpy as np  
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit
from scipy.integrate import quad

# -----------------------------
# Step 1: Load data
# -----------------------------
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v2.db"
)

run_id = 107
dataset = qc.load_by_id(run_id)
param_data = dataset.get_parameter_data()

# Manually specify keys
signal_key = "avg_avg_psd_nodrive"
freq_key = "freq_param"

# Extract and sort data
y_data = param_data[signal_key][signal_key].flatten()
x_data = param_data[signal_key][freq_key].flatten() / 1e6  # Convert Hz to MHz
sorted_indices = np.argsort(x_data)
x_data = x_data[sorted_indices]
y_data = y_data[sorted_indices]

# -----------------------------
# Step 2: Define Lorentzian
# -----------------------------
def lorentzian(x, x0, gamma, A, y0):
    return y0 + (A * gamma**2) / ((x - x0)**2 + gamma**2)

# -----------------------------
# Step 3: Fit over entire data
# -----------------------------
x_fit_data = x_data
y_fit_data = y_data

# Initial guess for fitting
x0_guess = x_fit_data[np.argmax(y_fit_data)]
gamma_guess = 0.00003
A_guess = max(y_fit_data)
y0_guess = min(y_fit_data)
p0 = [x0_guess, gamma_guess, A_guess, y0_guess]

# Fit on all data
popt, pcov = curve_fit(lorentzian, x_fit_data, y_fit_data, p0=p0)
perr = np.sqrt(np.diag(pcov))  # Standard errors

# Unpack fit results and uncertainties
x0_fit, gamma_fit, A_fit, y0_fit = popt
x0_err,  gamma_err,  A_err,  y0_err  = perr

print("\n✅ Fit parameters (±1σ uncertainties):")
print(f"Center frequency x0     = {x0_fit:.6f} ± {x0_err:.2e} MHz")
print(f"Width (gamma)           = {gamma_fit:.6f} ± {gamma_err:.2e} MHz")
print(f"Amplitude (A)           = {A_fit:.3e} ± {A_err:.1e} W/Hz")
print(f"Offset (y0)             = {y0_fit:.3e} ± {y0_err:.1e} W/Hz")

# -----------------------------
# Step 4: Area under the peak
# -----------------------------
lorentz_fit = lambda x: lorentzian(x, *popt)
background = lambda x: y0_fit

x1 = x_data.min()
x2 = x_data.max()

area_total, _ = quad(lorentz_fit, x1, x2)
area_background, _ = quad(background, x1, x2)
net_area = area_total - area_background

print(f"\n📐 Area under Lorentzian from {x1:.4f} to {x2:.4f} MHz:")
print(f"  Total area (with offset):     {area_total:.3e} W/Hz·MHz")
print(f"  Subtracted offset area:       {area_background:.3e} W/Hz·MHz")
print(f"  ✅ Net peak area:              {net_area:.3e} W/Hz·MHz")

# -----------------------------
# Step 4.5: Estimate error on area (Monte Carlo)
# -----------------------------
N_samples = 1000
rng = np.random.default_rng(seed=42)
samples = rng.multivariate_normal(popt, pcov, size=N_samples)
area_samples = []

for params in samples:
    area = quad(lambda x: lorentzian(x, *params), x1, x2)[0]
    offset = quad(lambda x: params[3], x1, x2)[0]
    area_samples.append(area - offset)

area_samples = np.array(area_samples)
area_mean = np.mean(area_samples)
area_std = np.std(area_samples)

print(f"\n📊 Monte Carlo estimate of net peak area:")
print(f"  Mean net area:         {area_mean:.3e} W/Hz·MHz")
print(f"  Std deviation (±1σ):   {area_std:.3e} W/Hz·MHz")

# -----------------------------
# Step 5: Plot with shaded area
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x_data, y_data, 'o', color='black', markersize=3, label="Raw PSD Data")

x_dense = np.linspace(x1, x2, 1000)
y_lorentz = lorentz_fit(x_dense)
y_offset = np.full_like(x_dense, y0_fit)

plt.plot(x_dense, y_lorentz, '-', color='red', linewidth=2, label="Lorentzian Fit (full range)")
plt.fill_between(x_dense, y_lorentz, color='yellow', alpha=0.4, label='Total area (with offset)')
plt.fill_between(x_dense, y_offset, y_lorentz, color='orange', alpha=0.6, label='Net peak area (offset subtracted)')

plt.xlabel("Frequency (MHz)", fontsize=14)
plt.ylabel(signal_key + " (W/Hz)", fontsize=14)
plt.title(f"Run {run_id} — Lorentzian Fit with Error and Area", fontsize=16)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

