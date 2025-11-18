import math
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
import numpy as np
from matplotlib.colors import LogNorm

# -----------------------
# User settings
# -----------------------
# Database location (optional)
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v30_14_11_25.db"
)

# Datasets to load (put all the IDs you want here)
#dataset_ids = [525,528,959, 956,864,897,971,900,906,909,912,918,940,968,962,959,956,925,953,946,971,928,915]  # e.g. [959, 956, 960, 961, ...]

dataset_ids = [
    525, 528, 864, 897, 900,
    906, 909, 912, 915, 918,
    925, 928, 940, 946, 953,
    956, 956, 959, 959, 962,
    965,968, 971, 971
]

# Names of parameters in your data
signal_name = "signal_shift_Vx_deriv_trc"
x_param_name = "QDAC_ch03_dc_constant_V"   # X-axis parameter
y_param_name = "QDAC_ch02_dc_constant_V"   # Y-axis parameter

# -----------------------
# Helper function
# -----------------------
def load_dataset_2d(ds_id, signal_name, x_name, y_name):
    """
    Load one dataset, return x_values, y_values, and G grid.

    Flips G vertically if y was scanned from high to low.
    """
    dataset = qc.load_by_id(ds_id)

    # Get conductance (or whatever signal you want)
    pdf_temp = dataset.to_pandas_dataframe_dict()
    G_raw = pdf_temp[signal_name]
    G_np = np.array(G_raw)

    # Get parameter data (we assume one main dependent parameter)
    interdeps = dataset.description.interdeps
    param_spec = interdeps.non_dependencies[0]
    param_name = param_spec.name

    data_xy = dataset.get_parameter_data(param_spec)

    x_raw = data_xy[param_name][x_name]
    y_raw = data_xy[param_name][y_name]

    x_np = np.array(x_raw)
    y_np = np.array(y_raw)

    # Unique x and y
    x_vals = np.unique(x_np)
    y_vals = np.unique(y_np)

    # Build 2D grid
    G_grid = np.zeros((len(y_vals), len(x_vals)))

    for m in range(len(y_vals)):
        for n in range(len(x_vals)):
            idx = m * len(x_vals) + n
            if idx < G_np.size:
                G_grid[m, n] = G_np[idx]

    # Check direction of y scan and flip if needed
    if len(y_np) > 1 and y_np[0] > y_np[-1]:
        print(f"Dataset {ds_id}: Y-axis scanned high→low. Flipping vertically.")
        G_grid = np.flipud(G_grid)

    # Sort y for plotting
    y_vals = np.sort(y_vals)

    return x_vals, y_vals, G_grid

# -----------------------
# Load all datasets
# -----------------------
x_list = []
y_list = []
G_list = []

for ds_id in dataset_ids:
    print(f"Loading dataset {ds_id}...")
    x_vals, y_vals, G_grid = load_dataset_2d(
        ds_id,
        signal_name=signal_name,
        x_name=x_param_name,
        y_name=y_param_name,
    )
    x_list.append(x_vals)
    y_list.append(y_vals)
    G_list.append(G_grid)

print("All datasets loaded.")

# -----------------------
# Plot individual datasets
# -----------------------
#n = len(dataset_ids)
#cols = math.ceil(math.sqrt(n))
#rows = math.ceil(n / cols)

#plt.figure(1, figsize=(4 * cols, 4 * rows))

#for i, ds_id in enumerate(dataset_ids):
#    plt.subplot(rows, cols, i + 1)
#    plt.pcolor(x_list[i], y_list[i], G_list[i])
#    plt.colorbar(label="Conductance G")
#    plt.xlabel(x_param_name)
#    plt.ylabel(y_param_name)
#    plt.title(f"Dataset {ds_id}")

#plt.tight_layout()
#plt.show()

#points of Marta's measurement October '25
g2_y_M = [0.389, 0.3732, 0.3575, 0.3417, 0.326, 0.3102]
g3_x_M = [0.318, 0.3363, 0.3543, 0.3723, 0.390, 0.408]

#points of mesaurement 181125
g2_y_S1 = [0.585, 0.535, 0.470, 0.4025, 0.310, 0.250]
g3_x_S1 = [0, 0.100, 0.200, 0.300, 0.400, 0.450]

# -----------------------
# Stitched plot (all datasets overlaid)
# -----------------------
# Global x/y limits
x_combined = np.concatenate(x_list)
y_combined = np.concatenate(y_list)
x_min, x_max = np.min(x_combined), np.max(x_combined)
y_min, y_max = np.min(y_combined), np.max(y_combined)

plt.figure(2, figsize=(8, 6))

# Try log-like plotting; fall back to linear if it fails

# If you truly want log scale, you can uncomment LogNorm and ensure data > 0
# norm = LogNorm(vmin=np.nanmax([np.min(G[G>0]) for G in G_list]), vmax=np.max([np.max(G) for G in G_list]))
for x_vals, y_vals, G_grid in zip(x_list, y_list, G_list):
        plt.pcolor(x_vals, y_vals, G_grid, alpha=0.6, cmap="plasma",
           vmin=0, vmax=100e-6)
  # , norm=norm)
plt.colorbar(label="peak shift estimate [V]")  # adjust label if using LogNorm

plt.scatter(g3_x_M, g2_y_M, marker='x', color='red', s=50)
plt.scatter(g3_x_S1, g2_y_S1, marker='x', color='green', s=50)
plt.xlabel(x_param_name)
plt.ylabel(y_param_name)
plt.title("Stitched 2D Plot: All datasets")
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.tight_layout()
plt.show()

print("\nProcessing complete!")
