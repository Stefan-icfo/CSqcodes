import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc

# Function to calculate moving average (not used in this code, but kept for potential use)
def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# Database location
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_zurichdata.db"
)

# Load experiments
experiments = qc.experiments()

# Load dataset
dataset_temp = qc.load_by_id(174)
df_temp = dataset_temp.to_pandas_dataframe_dict()

# Extract Hall resistance data
interdeps = dataset_temp.description.interdeps
param_spec = interdeps.non_dependencies[0]  # Hall resistance data
data_x = dataset_temp.get_parameter_data(param_spec)

# Get the trace data
trace = np.array(df_temp["v_r"])  # Make sure 'v_r' is a valid column in the DataFrame

# Ensure trace is valid
if trace.size == 0:
    raise ValueError("No data found in 'v_r'.")

num_points = len(trace)  # Get the number of points in v_r data
time_array = np.linspace(2, 5, num_points)  # Create time array from 2 to 5 seconds

# Save time_array and trace to a .txt file
output_filename = "time_trace_data.txt"
with open(output_filename, "w") as f:
    f.write("Time (s)\tTrace (µV)\n")  # Header
    for t, v in zip(time_array, trace.flatten()):  # Use flatten to ensure 1D
        f.write(f"{t:.6f}\t{v:.6f}\n")  # Write each pair of time and trace values

print(f"Data saved to {output_filename}")

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(time_array, trace, marker='o', linestyle='-', color='b', label='v_r data')
plt.title('Hall Resistance Data')
plt.xlabel('Time (s)')
plt.ylabel('v_r (µV)')
plt.grid(True)
plt.legend()
plt.show()

