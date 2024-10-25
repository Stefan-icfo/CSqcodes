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

# Step 2: Create a figure and axis for the main plot
fig, ax1 = plt.subplots(figsize=(10, 6))
line1, = ax1.plot(voltage_array, G_raw, marker='o', linestyle='-', color='b', label='G_raw vs Voltage')
ax1.set_title('Conductance (G_raw) vs Voltage')
ax1.set_xlabel('Voltage (V)')
ax1.set_ylabel('Conductance (S)')
ax1.grid()
ax1.legend()

# Variables to hold the selected x-values and the area under the curve
selected_x1 = None
selected_x2 = None
calculated_area = None

# Step 3: Function to calculate area under the selected region
def calculate_area(x1, x2):
    global selected_x1, selected_x2, calculated_area
    
    # Store the selected x-values
    selected_x1, selected_x2 = x1, x2
    
    # Find the indices of the selected x-values
    indices = np.where((voltage_array >= x1) & (voltage_array <= x2))
    selected_voltages = voltage_array[indices]
    selected_G = G_raw[indices]
    
    # Calculate area using the trapezoidal rule
    calculated_area = np.trapz(selected_G, selected_voltages)
    print(f'Area under the curve between {x1:.4f} V and {x2:.4f} V: {calculated_area:.4e} S·V')

    # Shade the area under the curve
    ax1.fill_between(selected_voltages, selected_G, color='yellow', alpha=0.5)

    # Redraw the figure to update the shaded area
    plt.draw()

# Step 4: Function to be called on selection
def on_select(eclick, erelease):
    x1, x2 = eclick.xdata, erelease.xdata
    calculate_area(x1, x2)
    plt.close(fig)  # Close the first plot after selection

# Step 5: Create RectangleSelector
selector = RectangleSelector(ax1, on_select, drawtype='box', useblit=True,
                             button=[1],  # Left mouse button
                             minspanx=5, minspany=5,
                             spancoords='pixels',
                             interactive=True)

# Step 6: Adjust Y-axis limits for the first plot
ax1.set_ylim(bottom=0, top=0.0001)  # Adjust this range as needed

plt.show()

# Step 7: Create a second plot after closing the first plot
if calculated_area is not None and selected_x1 is not None and selected_x2 is not None:
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    
    # Shade the area under the selected region in the second plot
    indices = np.where((voltage_array >= selected_x1) & (voltage_array <= selected_x2))
    selected_voltages = voltage_array[indices]
    selected_G = G_raw[indices]

    ax2.plot(voltage_array, G_raw, color='b', label='G_raw vs Voltage')
    ax2.fill_between(selected_voltages, selected_G, color='yellow', alpha=0.5)

    ax2.set_title('Area Under the Peak')
    ax2.set_xlabel('Voltage (V)')
    ax2.set_ylabel('Conductance (S)')
    ax2.grid()
    ax2.legend()
    ax2.set_ylim(bottom=0, top=0.0001)  # Adjust this range as needed

    # Display the calculated area on the plot
    ax2.text(0.5, 0.9, f'Area: {calculated_area:.4e} S·V', fontsize=12, ha='center', transform=ax2.transAxes)

    plt.show()
