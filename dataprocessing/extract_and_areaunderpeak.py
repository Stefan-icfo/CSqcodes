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

# Plot the raw data before asking for the area
plt.figure(figsize=(10, 6))
plt.plot(time_array_avg, trace_avg, marker='o', linestyle='-', color='b', label='Averaged v_r data (3-point)')
plt.title('Hall Resistance Data (Averaged)')
plt.xlabel('Time (s)')
plt.ylabel('v_r (µV)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Ask user for x1, x2 points and the y offset
run_id=1134 # Get the run_id from the description
print(f"\nSelect the area to measure for the Run ID {run_id}:")
x1 = float(input("Enter the x1 point (time in seconds): "))
x2 = float(input("Enter the x2 point (time in seconds): "))
y_offset = float(input("Enter the y offset value above which to calculate the area: "))

# Convert x1, x2 to indices
x1_idx = np.searchsorted(time_array_avg, x1)
x2_idx = np.searchsorted(time_array_avg, x2)

# Ensure valid indices
if x1_idx >= x2_idx:
    raise ValueError("x1 must be less than x2!")

# Calculate the area under the curve (integrate using the trapezoidal rule)
area = np.trapz(trace_avg[x1_idx:x2_idx] - y_offset, time_array_avg[x1_idx:x2_idx])

# Plot the data and highlight the area
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(time_array_avg, trace_avg, label='Data')
ax.fill_between(time_array_avg[x1_idx:x2_idx], trace_avg[x1_idx:x2_idx] - y_offset, color='yellow', alpha=0.6, label='Area Under Curve')

ax.set_xlabel('Time (s)', fontsize=14)
ax.set_ylabel('Phase Data (µV)', fontsize=14)
ax.legend()
ax.grid(True)
ax.set_title(f"Area under the curve for Run ID {run_id}")

plt.tight_layout()
plt.show()

# Print the calculated area
print(f"Calculated area under the curve from {x1} to {x2}, above y={y_offset}: {area:.6f} µV·s")








