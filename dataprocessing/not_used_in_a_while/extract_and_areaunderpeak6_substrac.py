import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit
from scipy.integrate import quad

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
x, y = load_sorted_data(989) #data
x_B, y_B = load_sorted_data(991) #beakground needed to be substracted

# Check if x-axes match
assert np.allclose(x, x_B), "X-axes do not match between the two runs!"

# -----------------------------
# Subtract spectra
# -----------------------------
y_diff = y - y_B

# -----------------------------
# Plot raw, background, and subtracted data
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x, y, color='black', label='Run 904 (raw)')
plt.plot(x_B, y_B, color='gold', label='Run 902 (background)')
plt.plot(x, y_diff, color='green', label='Subtracted (904 - 902)')
plt.xlabel("Frequency (MHz)", fontsize=14)
plt.ylabel("PSD (W/Hz)", fontsize=14)
plt.title("Spectrum with background subtracted", fontsize=16)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# -----------------------------
# Lorentzian fit on subtracted data
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
# Area under Lorentzian
# -----------------------------
lorentz_fit = lambda x: lorentzian(x, *popt)
background = lambda x: y0_fit

x1, x2 = x_fit.min(), x_fit.max()
area_total, _ = quad(lorentz_fit, x1, x2)
area_background, _ = quad(background, x1, x2)
net_area = area_total - area_background

print(f"\n📐 Area under Lorentzian:")
print(f"  Total area     = {area_total:.3e} W/Hz·MHz")
print(f"  Offset area    = {area_background:.3e} W/Hz·MHz")
print(f"  ✅ Net peak area = {net_area:.3e} W/Hz·MHz")

# -----------------------------
# Monte Carlo estimate of area uncertainty
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

print(f"\n📊 Monte Carlo estimate of area:")
print(f"  Mean area     = {area_mean:.3e} W/Hz·MHz")
print(f"  Std deviation = {area_std:.3e} W/Hz·MHz")

# -----------------------------
# Final plot with Lorentzian fit
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x, y_diff, 'o', color='green', markersize=3, label="Subtracted data")

x_dense = np.linspace(x1, x2, 1000)
y_lorentz = lorentz_fit(x_dense)
y_offset = np.full_like(x_dense, y0_fit)

plt.plot(x_dense, y_lorentz, '-', color='red', linewidth=2, label="Lorentzian fit")
plt.fill_between(x_dense, y_lorentz, color='lightgreen', alpha=0.4, label='Total area')
plt.fill_between(x_dense, y_offset, y_lorentz, color='lime', alpha=0.5, label='Net peak area')

plt.xlabel("Frequency (MHz)", fontsize=14)
plt.ylabel("PSD (W/Hz)", fontsize=14)
plt.title("Lorentzian fit on subtracted spectrum", fontsize=16)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

