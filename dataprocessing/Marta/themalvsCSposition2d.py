import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import re

# ---------------------
# CONFIGURATION
# ---------------------
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v3.db"
)

background_id = 322
run_ids = list(range(324, 372, 2))  # Even run IDs from 132 to 320

signal_key = "avg_avg_psd_nodrive"
freq_key = "freq_param"

# ---------------------
# LOAD PSD FUNCTION
# ---------------------
def load_psd(run_id):
    ds = qc.load_by_id(run_id)
    param_data = ds.get_parameter_data()
    block = param_data[signal_key]
    y = block[signal_key].flatten()
    x = block[freq_key].flatten() / 1e6  # Hz → MHz
    idx = np.argsort(x)
    return x[idx], y[idx], ds

# ---------------------
# PARSE GATE VOLTAGE
# ---------------------
def extract_gcs_voltage(name):
    match = re.search(r"gcs=([\d.]+)", name)
    return float(match.group(1)) if match else None

# ---------------------
# LOAD BACKGROUND
# ---------------------
x_bg, P_bg, _ = load_psd(background_id)
V_bg = np.sqrt(P_bg)

# ---------------------
# PROCESS MEASUREMENTS
# ---------------------
spectra = []
vgates = []
lengths = []

for run_id in run_ids:
    try:
        x, P, ds = load_psd(run_id)
        V = np.sqrt(P)
        V_interp_bg = np.interp(x, x_bg, V_bg)
        V_corr = V - V_interp_bg
        P_corr = V_corr**2

        vg = extract_gcs_voltage(ds.exp_name)
        if vg is None:
            print(f"⚠️ Skipping {run_id}: V_gcs not found in name '{ds.exp_name}'")
            continue

        spectra.append(P_corr)
        vgates.append(vg)
        lengths.append(len(P_corr))
        print(f"✔️ Loaded run {run_id} at V_gcs = {vg:.3f} mV")

    except Exception as e:
        print(f"❌ Error with run {run_id}: {e}")

# ---------------------
# PREPARE 2D GRID
# ---------------------
min_len = min(lengths)
freqs = x[:min_len]
spectra = np.array([s[:min_len] for s in spectra])
vgates = np.array(vgates)

# Sort by gate voltage
sort_idx = np.argsort(vgates)
spectra = spectra[sort_idx]
vgates = vgates[sort_idx]

# Create image extent
extent = [freqs[0], freqs[-1], vgates[0], vgates[-1]]

# ---------------------
# PLOT
# ---------------------
plt.figure(figsize=(10, 6))
plt.imshow(spectra, aspect='auto', extent=extent, origin='lower', cmap='viridis')
plt.colorbar(label='Peak PSD Intensity (W/Hz)')
plt.xlabel("Peak Frequency (MHz)")
plt.ylabel("Gate Voltage V_gcs (mV)")
plt.title("350udrive Peak Map (Background Subtracted)")
plt.tight_layout()
plt.show()







