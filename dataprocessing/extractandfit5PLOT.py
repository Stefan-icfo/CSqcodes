import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
from scipy.optimize import curve_fit
import os

# Function to calculate moving average of n points
def moving_average(a, n=3):
    return np.convolve(a, np.ones(n)/n, mode='valid')

# Lorentzian function for fitting
def lorentzian(x, x0, gamma, A, y0):
    return y0 + (A * gamma**2) / ((x - x0)**2 + gamma**2)

# Database location
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part2.db"
)

# Load experiments
experiments = qc.experiments()

# Load dataset
dataset_temp = qc.load_by_id(1134)
df_temp = dataset_temp.to_pandas_dataframe_dict()

# Extract Hall resistance data
interdeps = dataset_temp.description.interdeps
param_spec = interdeps.non_dependencies[0]  # Hall resistance data
data_x = dataset_temp.get_parameter_data(param_spec)

# Get the trace data
trace = np.array(df_temp["I_rf"])  # Make sure 'v_r' is a valid column in the DataFrame

# Ensure trace is valid
if trace.size == 0:
    raise ValueError("No data found in 'v_r'.")

num_points = len(trace)  # Get the number of points in v_r data
time_array = np.linspace(399, 402, num_points)  # Create time array from 2 to 5 seconds

# Calculate moving average with 3 points
trace_avg = moving_average(trace.flatten(), n=40)
time_array_avg = time_array[:len(trace_avg)]  # Adjust time array for averaged data

# Save time_array and trace_avg to a .txt file
output_filename = r"C:\\Users\\LAB-nanooptomechanic\\Desktop\\ringdown\\frequency.txt"
with open(output_filename, "w") as f:
    f.write("Time (s)\tTrace (µV)\n")  # Header
    for t, v in zip(time_array_avg, trace_avg):
        f.write(f"{t:.6f}\t{v:.6f}\n")  # Write each pair of time and trace values

print(f"Data saved to {output_filename}")

# Fit the data to a Lorentzian curve
p0 = [4.8, 0.01, max(trace_avg), min(trace_avg)]  # Initial guess for [x0, gamma, A, y0]
popt, pcov = curve_fit(lorentzian, time_array_avg, trace_avg, p0=p0)

# Generate fitted Lorentzian curve
time_fit = np.linspace(min(time_array_avg), max(time_array_avg), 1000)
trace_fit = lorentzian(time_fit, *popt)

# Plot the averaged data and the Lorentzian fit
plt.figure(figsize=(10, 6))
plt.plot(time_array_avg, trace_avg, marker='o', linestyle='-', color='b', label='Averaged v_r data (3-point)')
#plt.plot(time_fit, trace_fit, linestyle='--', color='r', label='Lorentzian Fit')
plt.title('Hall Resistance Data (Averaged with Lorentzian Fit)')
plt.xlabel('Time (s)')
plt.ylabel('v_r (µV)')

plt.show()

# Print fitted parameters
print("Fitted parameters:")
print(f"x0 (center): {popt[0]:.6f}")
print(f"gamma (width): {popt[1]:.6f}")
print(f"A (amplitude): {popt[2]:.6f}")
print(f"y0 (offset): {popt[3]:.6f}")

