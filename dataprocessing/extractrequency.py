import qcodes as qc
import matplotlib.pyplot as plt
import numpy as np
import os

# -------------------- QCoDeS Setup --------------------
qc.config["core"]["db_location"] = (
    r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD11_D7_C1_part2.db"
)

# -------------------- Load Dataset --------------------
run_id = 910
ds = qc.load_by_id(run_id)
df = ds.to_pandas_dataframe()

# -------------------- Extract Data --------------------
# Convert zurich_osc0_freq from Hz to MHz
x = df["zurich_osc0_freq"].values / 1e6  # x-axis: frequency in MHz
y = df["I_rf_avg"].values                # y-axis: I_rf_avg signal

# -------------------- Plot Raw Data --------------------
plt.figure(figsize=(10, 6))
plt.plot(x, y, 'o', markersize=4, label='Raw I_rf_avg')
plt.xlabel("Frequency (MHz)")
plt.ylabel("I_rf_avg (µV)")
plt.title(f"Run ID: {ds.run_id}")
plt.grid(True)
plt.legend()
plt.show()

# -------------------- 3-Point Averaging --------------------
# Trim to a multiple of 3
remainder = len(x) % 3
if remainder:
    x = x[:-remainder]
    y = y[:-remainder]

# Compute averages
x_avg = x.reshape(-1, 3).mean(axis=1)
y_avg = y.reshape(-1, 3).mean(axis=1)

# -------------------- Plot Averaged Data --------------------
plt.figure(figsize=(10, 6))
plt.plot(x, y, 'o', color='gray', alpha=0.3, label='Raw I_rf_avg')
plt.plot(x_avg, y_avg, 'o-', color='orange', label='3-Point Average')
plt.xlabel("Frequency (MHz)")
plt.ylabel("I_rf_avg (µV)")
plt.title(f"3-Point Averaged I_rf_avg | Run ID: {ds.run_id}")
plt.grid(True)
plt.legend()
plt.show()

# -------------------- Save Averaged Data --------------------
output_path = r"C:\Users\LAB-nanooptomechanic\Desktop\run910_Irfavg.txt"
os.makedirs(os.path.dirname(output_path), exist_ok=True)

with open(output_path, 'w') as f:
    f.write("Frequency (MHz)\tI_rf_avg (µV)\n")
    for xi, yi in zip(x_avg, y_avg):
        f.write(f"{xi:.6f}\t{yi:.6f}\n")

print(f"Averaged data saved to: {output_path}")





