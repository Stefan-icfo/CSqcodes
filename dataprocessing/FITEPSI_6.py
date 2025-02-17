import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import Boltzmann, Planck, h, e
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
from scipy.optimize import curve_fit

# Constants
kB = Boltzmann  # Boltzmann constant in Hz/K (for angular frequency units)
Te_fixed = 0.2 # Fixed temperature in Kelvin
t_initial =6e9 * h  # Initial guess for t in Hz
alpha_fixed = 0.37  # Fixed alpha value in eV/V
xo_initial = 0.01  # Initial guess for xo in the same units as delta
qc.config["core"]["db_location"] = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part3"
# Bounds for fitting parameters


# File numbers for loading the datafile_numbers = [num for num in range(3011,3018) if num not in [3016,3015,3017] ]
file_numbers = [646]
# Lists to store t for all measurements
t_values_GHz = []

# Normalization function
def normalize(M):
    M_max = np.mean(M[:20])  # Average of the first points for M = 1
    M_min = np.mean(M[-30:])  # Average of the last points for M = 0
    M_normalized = (M - M_min) / (M_max - M_min)
    return M_normalized

# Function to find the delta shift where M is closest to 0.5
def center_curve(delta, M_norm):
    idx_closest = np.argmin(np.abs(M_norm - 0.5))
    delta_shifted = delta - delta[idx_closest]
    return delta_shifted

# Define the function f_epsilon_slope_fit with Te as a fixed parameter, alpha as fixed, and xo as a free parameter
def f_epsilon_slope_fit(average_delta_shifted, t, xo):
    term = (average_delta_shifted - xo) * alpha_fixed * e
    return (1/2) * (1 - (term / np.sqrt(term**2 + (4 * t**2))) * np.tanh(np.sqrt(term**2 + (4 * t**2)) / (2 * kB * Te_fixed)))

# Loop over each measurement file
for num in file_numbers:
    # Load data
    delta = np.load(f'meas{num}_delta_array.npy')
    peakpos = np.load(f'meas{num}_peakpos_V.npy')

    # Normalize and shift data
    M_norm = normalize(peakpos)
    delta_shifted = center_curve(delta, M_norm)

    # Perform the curve fitting
    popt, pcov = curve_fit(f_epsilon_slope_fit, delta_shifted, M_norm, p0=[t_initial, xo_initial])

    # Extract the fitted parameters and their standard errors
    t_fit, xo_fit = popt
    t_error, xo_error = np.sqrt(np.diag(pcov))

    # Convert t to GHz
    t_fit_GHz = t_fit / Planck / 1e9
    t_error_GHz = t_error / Planck / 1e9

    # Store t (in GHz) for averaging later
    t_values_GHz.append(t_fit_GHz)

    # Plot the data and the fit
    plt.figure(figsize=(10, 6))
    plt.plot(delta_shifted, M_norm, 'b+', label=f'Data for meas{num}')
    plt.plot(delta_shifted, f_epsilon_slope_fit(delta_shifted, *popt), 'r-', label=f'Fit: t={t_fit_GHz:.3f} ± {t_error_GHz:.3f} GHz, xo={xo_fit:.3e} ± {xo_error:.3e}')
    plt.xlabel('Shifted Delta (V)')
    plt.ylabel('Normalized M')
    plt.title(f'Fit of f_epsilon_slope_fit Model to Data (Measurement {num})')
    plt.legend()
    plt.show()

    # Print the results for each measurement
    print(f"Measurement {num}:")
    print(f"Fitted t: {t_fit_GHz:.3f} GHz ± {t_error_GHz:.3f} GHz")
    print(f"Fitted xo: {xo_fit:.3e} ± {xo_error:.3e}")
    print("----------------------------------------------------")

# Calculate average and standard deviation for t (in GHz)
t_avg_GHz = np.mean(t_values_GHz)
t_std_GHz = np.std(t_values_GHz)

# Print the average and standard deviation of t
print("Summary of Fitted Parameters:")
print(f"Average t (in GHz): {t_avg_GHz:.3f} GHz ± {t_std_GHz:.3f} GHz")
