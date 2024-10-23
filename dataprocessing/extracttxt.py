import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# Set the database location
qc.config["core"]["db_location"] = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1.db"

# Load the experiment
experiment_id = 314
dataset = qc.load_by_id(experiment_id)
df = dataset.to_pandas_dataframe_dict()

# Extract data
trace = np.array(df["G"])

# Smooth the data
#smooth_data = moving_average(trace, 30)

# Generate frequency axis
voltage = np.linspace(-1, -0.3, len(trace))

# Plot
plt.plot(voltage, trace)
plt.xlabel("voltage (Hz)")
plt.ylabel("G (a.u.)")
plt.show()

# Save the data to a text file
output_file_path = "Z:\\Users\\Marta\\gvg\\314.txt"
with open(output_file_path, 'w') as f:
    for freq, amp in zip(voltage, trace):
        f.write(f"{voltage}, {trace}\n")

print(f"Data extracted and saved to '{output_file_path}'")