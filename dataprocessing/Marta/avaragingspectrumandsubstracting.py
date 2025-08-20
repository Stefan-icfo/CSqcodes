import numpy as np 
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit

# -----------------------------
# CONFIGURATION
# -----------------------------
run_data = 162
run_background = 165

qc.config["core"]["db_location"] = (
    "C:/Users/LAB-nanooptomechanic/Documents/MartaStefan/CSqcodes/Data/Raw_data/CD12_B5_F4v2.db"
)

# -----------------------------
# Lorentzian model with y0 = 0
# -----------------------------
def lorentzian_fixed_y0(x, x0, gamma, A):
    return (A * gamma**2) / ((x - x0)**2 + gamma**2)

# -----------------------------
# Load PSD from run_data
# -----------------------------
dataset = qc.load_by_id(run_data)
data = dataset.get_parameter_data('avg_psd')['avg_psd']
freq = np.array(data['freq_param'])  # in Hz
time = np.array(data['time_param'])
psd_raw = np.array(data['avg_psd'])

unique_time = np.unique(time)
unique_freq = np.unique(freq)
n_time = len(unique_time)
n_freq = len(unique_freq)

psd = psd_raw.reshape(n_time, n_freq)
freq = unique_freq
time = unique_time

# -----------------------------
# Identify good and bad traces
# -----------------------------
bad_mask = (time >= 4440) & (time <= 5000)
good_mask = ~bad_mask
psd_good = psd[good_mask]

# -----------------------------
# Plot 1: All traces (black/red)
# -----------------------------
plt.figure(figsize=(10, 6))
for i in range(n_time):
    color = 'red' if bad_mask[i] else 'black'
    plt.plot(freq, psd[i], color=color, alpha=0.07, linewidth=0.5)
plt.xlabel("Frequency (Hz)")
plt.ylabel("PSD (W/Hz)")
plt.title("All traces (black) and excluded (red)")
plt.grid(True)
plt.tight_layout()
plt.show()

# -----------------------------
# Plot 2: Mean with and without bad traces
# -----------------------------
mean_all = np.mean(psd, axis=0)
mean_good = np.mean(psd_good, axis=0)

plt.figure(figsize=(10, 6))
plt.plot(freq, mean_good, color='black', label="Mean (good only)")
plt.plot(freq, mean_all, color='red', linestyle='--', label="Mean (all)")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Averaged PSD (W/Hz)")
plt.title("Mean PSD: good vs all traces")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# -----------------------------
# Plot 3: Mean of good traces only
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(freq, mean_good, color='black')
plt.xlabel("Frequency (Hz)")
plt.ylabel("Averaged PSD (W/Hz)")
plt.title("Mean PSD (good traces only)")
plt.grid(True)
plt.tight_layout()
plt.show()

# -----------------------------
# Load background PSD (run_background)
# -----------------------------
bkg_data = qc.load_by_id(run_background).get_parameter_data()
x_bkg = bkg_data['avg_avg_psd_nodrive']['freq_param'].flatten()
y_bkg = bkg_data['avg_avg_psd_nodrive']['avg_avg_psd_nodrive'].flatten()
P_noise = np.interp(freq, x_bkg, y_bkg)

# -----------------------------
# Voltage-domain subtraction
# -----------------------------
V_total = np.sqrt(np.clip(mean_good, 0, None))
V_noise = np.sqrt(np.clip(P_noise, 0, None))
V_signal = V_total - V_noise
P_signal = V_signal**2  # Final corrected spectrum

# -----------------------------
# Plot 4: Good mean vs subtracted
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(freq, mean_good, color='black', label="Mean PSD (good)")
plt.plot(freq, P_signal, color='green', label="Mean - Background", linewidth=1.5)
plt.xlabel("Frequency (Hz)")
plt.ylabel("PSD (W/Hz)")
plt.title("Corrected mean PSD (green) vs original (black)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# -----------------------------
# Fit Lorentzian to corrected data (y0 = 0)
# -----------------------------
fit_successful = False
try:
    x0_guess = freq[np.argmax(P_signal)]
    gamma_guess = 3e3
    A_guess = max(P_signal)
    p0 = [x0_guess, gamma_guess, A_guess]
    popt, pcov = curve_fit(lorentzian_fixed_y0, freq, P_signal, p0=p0, maxfev=10000)
    fit_successful = True
    x0, gamma, A = popt
    Q = x0 / (2 * gamma)
except Exception as e:
    print(f"⚠️ Lorentzian fit failed: {e}")
    fit_successful = False

# -----------------------------
# Plot 5: Corrected PSD + Lorentzian fit + results
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(freq, P_signal, color='green', marker='+', linestyle='None', label="Corrected PSD (subtracted)")

if fit_successful:
    freq_dense = np.linspace(freq.min(), freq.max(), 1000)
    plt.plot(freq_dense, lorentzian_fixed_y0(freq_dense, *popt), 'r--', label='Lorentzian fit (y₀=0)')

    fit_text = (
        f"Fit results:\n"
        f"$x_0$ = {x0/1e6:.6f} MHz\n"
        f"$\gamma$ = {gamma:.0f} Hz\n"
        f"$Q$ = {Q:.1f}\n"
        f"$A$ = {A:.2e} W/Hz\n"
        f"$y_0$ = 0 (fixed)"
    )
    plt.text(0.02, 0.95, fit_text, transform=plt.gca().transAxes,
             verticalalignment='top', fontsize=10,
             bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

plt.xlabel("Frequency (Hz)")
plt.ylabel("PSD (W/Hz)")
plt.title("Final PSD (corrected) with Lorentzian fit (y₀ = 0)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
