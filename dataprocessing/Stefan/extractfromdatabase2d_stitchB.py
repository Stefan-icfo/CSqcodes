import math
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
import numpy as np
from matplotlib.colors import LogNorm

# Database location
qc.config["core"]["db_location"] = "C:\\" + "Users\\" + "LAB-nanooptomechanic\\" + "Documents\\" + "MartaStefan\\" + "CSqcodes\\" + "Data\\" + "Raw_data\\" + "CD12_B5_F4v11.db"

experiments = qc.experiments()

# Load datasets
dataset1 = qc.load_by_id(1898)
dataset2 = qc.load_by_id(1902)

# Convert to pandas dataframes
pdf_temp1 = dataset1.to_pandas_dataframe_dict()
pdf_temp2 = dataset2.to_pandas_dataframe_dict()

# Get conductance data
G1_raw = pdf_temp1['G']
G2_raw = pdf_temp2['G']
G1_np = np.array(G1_raw)
G2_np = np.array(G2_raw)

# Get parameter data
interdeps = dataset1.description.interdeps
param_spec = interdeps.non_dependencies[0]
param_name = param_spec.name

data_xy1 = dataset1.get_parameter_data(param_spec)
data_xy2 = dataset2.get_parameter_data(param_spec)

xy1 = data_xy1[param_name][param_name]
xy2 = data_xy2[param_name][param_name]

# Extract x and y coordinates for both datasets
x1_raw = data_xy1[param_name]['QDAC_ch06_dc_constant_V']
y1_raw = data_xy1[param_name]['QDAC_ch02_dc_constant_V']
x2_raw = data_xy2[param_name]['QDAC_ch06_dc_constant_V']
y2_raw = data_xy2[param_name]['QDAC_ch02_dc_constant_V']

# Convert to numpy arrays
x1_np = np.array(x1_raw)
y1_np = np.array(y1_raw)
x2_np = np.array(x2_raw)
y2_np = np.array(y2_raw)

# Get unique values for grid
x1 = np.unique(x1_np)
y1 = np.unique(y1_np)
x2 = np.unique(x2_np)
y2 = np.unique(y2_np)

# Check scan direction for y-axis
print(f"Dataset 1 - First few y values: {y1_np[:10]}")
print(f"Dataset 2 - First few y values: {y2_np[:10]}")

# Create 2D arrays
G1 = np.zeros([len(y1), len(x1)])
G2 = np.zeros([len(y2), len(x2)])

# Fill the arrays
for m in range(len(y1)):
    for n in range(len(x1)):
        G1[m, n] = G1_np[m * len(x1) + n]

for m in range(len(y2)):
    for n in range(len(x2)):
        G2[m, n] = G2_np[m * len(x2) + n]

# Check if y1 was scanned from high to low and flip if necessary
if len(y1_np) > 1 and y1_np[0] > y1_np[-1]:
    print("Dataset 1: Y-axis was scanned from high to low. Flipping vertically...")
    G1 = np.flipud(G1)  # Flip vertically (up-down)

# Sort y1 to ensure it's in ascending order for plotting
y1 = np.sort(y1)
# y2 should already be sorted from np.unique
y2 = np.sort(y2)

# Plot individual datasets
plt.figure(1, figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.pcolor(x1, y1, G1)
plt.colorbar(label='Conductance G')
plt.xlabel('QDAC_ch06_dc_constant_V')
plt.ylabel('QDAC_ch02_dc_constant_V')
plt.title('Dataset 1 (Corrected)')

plt.subplot(1, 2, 2)
plt.pcolor(x2, y2, G2)
plt.colorbar(label='Conductance G')
plt.xlabel('QDAC_ch06_dc_constant_V')
plt.ylabel('QDAC_ch02_dc_constant_V')
plt.title('Dataset 2')

plt.tight_layout()
plt.show()

# Stitching now
x_combined = np.concatenate([x1, x2])
y_combined = np.concatenate([y1, y2])
x_min, x_max = np.min(x_combined), np.max(x_combined)
y_min, y_max = np.min(y_combined), np.max(y_combined)

# Plot both datasets with transparency and log scale
plt.figure(2, figsize=(8, 6))

try:
    plt.pcolor(x1, y1, G1, alpha=0.7, cmap='plasma', norm=LogNorm())
    plt.pcolor(x2, y2, G2, alpha=0.7, cmap='plasma', norm=LogNorm())
    plt.colorbar(label='Conductance G (log scale)')
except:
    print("Cannot use log scale (possibly due to negative/zero values), using linear scale")
    plt.pcolor(x1, y1, G1, alpha=0.7, cmap='plasma')
    plt.pcolor(x2, y2, G2, alpha=0.7, cmap='plasma')
    plt.colorbar(label='Conductance G')

plt.xlabel('QDAC_ch06_dc_constant_V')
plt.ylabel('QDAC_ch02_dc_constant_V')
plt.title('Stitched 2D Plot: Dataset 1 + Dataset 2 (Corrected)')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.show()

print("\nProcessing complete!")