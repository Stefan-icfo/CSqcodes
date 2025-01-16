import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
from utils.CS_utils import centered_moving_average
from scipy.optimize import curve_fit

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def lorentzian(f, f0, gamma, A, offset):
    return A * (gamma**2 / ((f - f0)**2 + gamma**2)) + offset

# Database location
qc.config["core"]["db_location"] = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part2.db"

experiments = qc.experiments()

dataset_temp = qc.load_by_id(1134)
df_temp = dataset_temp.to_pandas_dataframe_dict()

interdeps = dataset_temp.description.interdeps
param_spec = interdeps.non_dependencies[0]  # hall resistance data
data_x = dataset_temp.get_parameter_data(param_spec)

trace = np.array(df_temp["I_rf"])
num_points = len(trace)

# Frequency array from 399 MHz to 402 MHz
freq_array = np.linspace(399, 402, num_points)

# Compute moving average every 5 points
window_size = 30
averaged_trace = moving_average(trace, n=window_size)
averaged_freq = moving_average(freq_array, n=window_size)

# Initial guess for Lorentzian fit [f0, gamma, A, offset]
initial_guess = [averaged_freq[np.argmax(averaged_trace)], 0.01, max(averaged_trace), min(averaged_trace)]

# Fit the Lorentzian curve
popt, pcov = curve_fit(lorentzian, averaged_freq, averaged_trace, p0=initial_guess)

# Plot data and Lorentzian fit
plt.plot(averaged_freq, averaged_trace, color='blue')
plt.plot(averaged_freq, lorentzian(averaged_freq, *popt), 'r--', label='Lorentzian Fit')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Signal Intensity')
plt.title('Averaged Signal vs Frequency with Lorentzian Fit')
plt.show()
