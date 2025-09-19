import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from matplotlib.cm import get_cmap
from database import *

# -----------------------------
# Configuration
# -----------------------------
db_path = DATABASE_LOCATION #"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v2.db"
qc.config["core"]["db_location"] = db_path

start_id = 1323
end_id = 1371
default_signal_key = "avg_avg_psd_nodrive"
freq_key = "freq_param"

# -----------------------------
# Helper: Load PSD data safely
# -----------------------------
def load_psd(run_id):
    ds = qc.load_by_id(run_id)
    param_data = ds.get_parameter_data()

    if default_signal_key in param_data:
        key = default_signal_key
    else:
        key = list(param_data.keys())[0]
        print(f"⚠️ Run {run_id}: default key not found, using '{key}'.")

    try:
        y = param_data[key][key].flatten()
        x = param_data[key][freq_key].flatten() / 1e6
    except KeyError as e:
        raise KeyError(f"Missing subkeys in run {run_id}: {e}")

    if len(x) == 0 or len(y) == 0:
        raise ValueError(f"Run {run_id} has empty x or y data.")

    sorted_idx = np.argsort(x)
    return x[sorted_idx], y[sorted_idx], ds

# -----------------------------
# Collect all valid pairs
# -----------------------------
all_pairs = list(range(start_id, end_id, 4))
valid_traces = []

for i, bg_id in enumerate(all_pairs):
    data_id = bg_id + 2
    try:
        _, y_bg, _ = load_psd(bg_id)
        x_data, y_data, ds_data = load_psd(data_id)
        y_corr = y_data - y_bg

        # Extract Vg
        name = ds_data.exp_name
        if "gcs=" in name:
            try:
                v_gcs = float(name.split("gcs=")[1].split("mV")[0])
            except:
                v_gcs = data_id
        else:
            v_gcs = data_id

        valid_traces.append((x_data, y_corr, v_gcs))

    except Exception as e:
        print(f"❌ Skipping pair {bg_id}/{data_id}: {e}")
        continue

# -----------------------------
# Final overlay plot
# -----------------------------
if len(valid_traces) > 6:
    middle_traces = valid_traces[3:-3]

    # Normalize colormap from blue → red
    cmap = get_cmap("plasma")
    n = len(middle_traces)

    plt.figure(figsize=(10, 6))
    for j, (x, y, v_gcs) in enumerate(middle_traces):
        color = cmap(j / max(n - 1, 1))
        plt.plot(x, y, label=f"{v_gcs:.1f} mV", color=color)

    plt.xlabel("Frequency (MHz)", fontsize=14)
    plt.ylabel("PSD minus background (W/Hz)", fontsize=14)
    plt.title("Overlay (excluding first and last 3)", fontsize=16)
    plt.legend(title="Vg (mV)", fontsize=9)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
else:
    print("⚠️ Not enough valid traces for overlay (need > 6).")



