import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc

# Database location
qc.config["core"]["db_location"] = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part3.db"

# Load dataset
dataset = qc.load_by_id(391)

# Debugging: Check available parameters
print("Dataset parameters:", dataset.parameters)

# Fetch parameter data
data_dict = dataset.get_parameter_data()

# Debug: Print available keys in data_dict
print("Available parameter data keys:", data_dict.keys())

# Check if 'v_r' exists
if 'v_r' in data_dict:
    v_r = data_dict['v_r']['v_r']

    # Convert to NumPy array
    v_r = np.array(v_r).flatten()

    # Step 1: Get the number of points in v_r
    num_points = len(v_r)

    # Step 2: Create a new time array from 2 to 5 seconds with the same number of points
    time_range = np.linspace(2, 5, num_points)

    # Step 3: Plot v_r vs. time_range
    plt.figure(figsize=(10, 6))
    plt.plot(time_range, v_r, marker='o', linestyle='-', color='b', label="v_r vs. time_range")
    plt.xlabel("Time (s)")  # Label for X-axis
    plt.ylabel("v_r (V)")   # Label for Y-axis
    plt.title("v_r vs. Time (2 to 5 seconds)")
    plt.grid(True)
    plt.legend()
    plt.show()

else:
    print("Error: 'v_r' not found. Available keys:", data_dict.keys())









