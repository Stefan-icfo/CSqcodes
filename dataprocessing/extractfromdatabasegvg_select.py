import matplotlib.pyplot as plt 
import qcodes as qc
import numpy as np
import pandas as pd
from matplotlib.widgets import RectangleSelector

# Database location
qc.config["core"]["db_location"] = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1.db"
experiments = qc.experiments()

# Load the dataset by its ID
dataset_id = 2052  # Update this to your actual dataset ID
dataset = qc.load_by_id(dataset_id)

# Convert dataset to pandas DataFrame
pdf_2052 = dataset.to_pandas_dataframe_dict()
G_raw = pdf_2052["G"]  # Extract the G values

# Step 1: Create a new voltage array
start_vg = -1.234
stop_vg = -1.288
num_points = len(G_raw)  # Use the same number of points as G_raw

# Create an array of voltages linearly spaced
voltage_array = np.linspace(start_vg, stop_vg, num_points)

# Step 2: Plot G_raw against the new voltage array
fig, ax = plt.subplots(figsize=(10, 6))
line, = ax.plot(voltage_array, G_raw, marker='o', linestyle='-', color='b', label='G_raw vs Voltage')
plt.title('Conductance (G_raw) vs Voltage')
plt.xlabel('Voltage (V)')
plt.ylabel('Conductance (S)')
plt.grid()
plt.legend()

# Step 3: Function to calculate area under the selected region
def calculate_area(x1, x2):
    # Find the indices of the selected x-values
    indices = np.where((voltage_array >= x1) & (voltage_array <= x2))
    selected_voltages = voltage_array[indices]
    selected_G = G_raw[indices]
    
    # Calculate area using the trapezoidal rule
    area = np.trapz(selected_G, selected_voltages)
    print(f'Area under the curve between {x1:.4f} V and {x2:.4f} V: {area:.4e} S·V')

# Step 4: Function to be called on selection
def on_select(eclick, erelease):
    x1, x2 = eclick.xdata, erelease.xdata
    calculate_area(x1, x2)

# Step 5: Create RectangleSelector
selector = RectangleSelector(ax, on_select, drawtype='box', useblit=True,
                             button=[1],  # Left mouse button
                             minspanx=5, minspany=5,
                             spancoords='pixels',
                             interactive=True)

plt.show()
