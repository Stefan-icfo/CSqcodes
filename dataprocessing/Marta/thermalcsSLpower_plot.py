import numpy as np
import matplotlib.pyplot as plt

# Simulated frequency axis
freq = np.linspace(137.9, 138.1, 1000)  # MHz

# Simulated signal (Lorentzian) + noise floor
P_noise = 1e-7 * np.ones_like(freq)  # Flat background
P_total = P_noise + 1e-6 / (1 + ((freq - 138.0) / 0.005)**2)  # Lorentzian + noise

# Subtraction Method 1: Direct power subtraction
P_signal_direct = P_total - P_noise

# Subtraction Method 2: Voltage-domain subtraction
P_signal_voltage = (np.sqrt(P_total) - np.sqrt(P_noise))**2

# --- Plotting ---
plt.figure(figsize=(10, 6))
plt.plot(freq, P_total, 'k--', label="Total (signal + noise)", linewidth=1)
plt.plot(freq, P_noise, 'gray', linestyle=':', label="Noise background", linewidth=1)
plt.plot(freq, P_signal_direct, 'b-', label="Direct Subtraction", linewidth=2)
plt.plot(freq, P_signal_voltage, 'g-', label="Voltage-domain Subtraction", linewidth=2)

plt.title("Signal Extraction from PSD with Background Subtraction", fontsize=14)
plt.xlabel("Frequency (MHz)", fontsize=12)
plt.ylabel("PSD (W/Hz)", fontsize=12)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()




