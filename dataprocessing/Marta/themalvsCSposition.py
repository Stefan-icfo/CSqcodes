import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from matplotlib.cm import get_cmap

# -----------------------------
# Configuration
# -----------------------------
db_path = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v2.db"
qc.config["core"]["db_location"] = db_path

background_id = 469
measurement_ids = list(range(499, 500, 2))  # odd numbers only
signal_key = "avg_avg_psd_nodrive"
freq_key = "freq_param"

# -----------------------------
# Helper: Load data
# -----------------------------
def load_psd(run_id):
    ds = qc.load_by_id(run_id)
    param_data = ds.get_parameter_data()
    y = param_data[signal_key][signal_key].flatten()
    x = param_data[signal_key][freq_key].flatten() / 1e6  # Hz → MHz
    sorted_idx = np.argsort(x)
    return x[sorted_idx], y[sorted_idx], ds

# -----------------------------
# Step 1: Load background
# -----------------------------
x_bg, y_bg, _ = load_psd(background_id)
assert np.allclose(x_bg, x_bg), "Frequencies do not match between datasets!"

# -----------------------------
# Step 2: Loop through data runs
# -----------------------------
plt.figure(figsize=(10, 6))
cmap = get_cmap('tab10')
legends = []

for i, run_id in enumerate(measurement_ids):
    x, y, ds = load_psd(run_id)
    y_corrected = y - y_bg

    # Get gate voltage from name string, e.g. "thermalV_gcs=965 mV"
    name = ds.exp_name
    if "gcs=" in name:
        try:
            v_gcs = float(name.split("gcs=")[1].split("mV")[0])
        except Exception:
            v_gcs = run_id
    else:
        v_gcs = run_id

    plt.plot(x, y_corrected, label=f'{v_gcs:.1f} mV', color=cmap(i % 10))
    legends.append(f'{v_gcs:.1f} mV')

# -----------------------------
# Finalize plot
# -----------------------------
plt.xlabel("Frequency (MHz)", fontsize=14)
plt.ylabel("Corrected PSD (W/Hz)", fontsize=14)
plt.title("PSD vs Frequency (background subtracted)", fontsize=16)
plt.legend(title="Vg (gate voltage)", fontsize=10)
plt.grid(True)
plt.tight_layout()
plt.show()
