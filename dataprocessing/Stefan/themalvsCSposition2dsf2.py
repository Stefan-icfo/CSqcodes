#this code doesnt sort the spectra any longer, since with the vgcs adjustment it´s mixing things up - 211125

import os
import re
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ---------------------
# USER CONFIG
# ---------------------
qc.config["core"]["db_location"] = (
    r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v32_19_11_25.db"
)

background_id = 175
run_ids = list(range(177, 280, 3))  # change as needed

background_id = 293
run_ids = list(range(295, 398, 3))  # change as needed


signal_key = "avg_avg_psd_nodrive"
freq_key = "freq_param"

# ---------------------
# HELPER FUNCTIONS
# ---------------------
def load_psd(run_id):
    ds = qc.load_by_id(run_id)
    block = ds.get_parameter_data()[signal_key]
    y = block[signal_key].flatten()
    x = block[freq_key].flatten() / 1e6  # Hz → MHz
    idx = np.argsort(x)
    return x[idx], y[idx], ds

def extract_gcs_voltage(name):
    m = re.search(r"gcs_=([\d.]+)\s*mV", name)
    return float(m.group(1)) if m else None

# ---------------------
# LOAD BACKGROUND
# ---------------------
x_bg, P_bg, _ = load_psd(background_id)
V_bg = np.sqrt(P_bg)

# ---------------------
# PROCESS MEASUREMENTS (IN DB ORDER)
# ---------------------
spectra = []
vgates = []
lengths = []
valid_run_ids = []

for run_id in run_ids:
    try:
        x, P, ds = load_psd(run_id)
        V = np.sqrt(P)
        V_interp_bg = np.interp(x, x_bg, V_bg)
        V_corr = V - V_interp_bg
        P_corr = V_corr**2

        vg = extract_gcs_voltage(ds.exp_name)
        if vg is None:
            print(f"⚠️ Skipping {run_id}: V_gcs not found in '{ds.exp_name}'")
            continue

        spectra.append(P_corr)
        vgates.append(vg)
        lengths.append(len(P_corr))
        valid_run_ids.append(run_id)
        print(f"✔️ Loaded run {run_id} at V_gcs = {vg:.3f} mV")

    except Exception as e:
        print(f"❌ Error with run {run_id}: {e}")

if not spectra:
    raise RuntimeError("No valid spectra loaded – check run_ids / signal_key / regex.")

# ---------------------
# PREPARE 2D GRID (NO SORTING BY GATE)
# ---------------------
min_len = min(lengths)
freqs = x[:min_len]
spectra = np.array([s[:min_len] for s in spectra])
vgates = np.array(vgates)
valid_run_ids = np.array(valid_run_ids)

# x: frequency in MHz, y: trace index (0, 1, 2, ...)
extent = [freqs[0], freqs[-1], 0, spectra.shape[0] - 1]

# ---------------------
# PLOT (ORDER = DB ORDER)
# ---------------------
plt.figure(figsize=(10, 6))
plt.imshow(
    spectra,
    aspect="auto",
    extent=extent,
    origin="lower",
    cmap="viridis",
    interpolation="none",
)
plt.colorbar(label="Peak PSD Intensity (W/Hz)")
plt.xlabel("Peak Frequency (MHz)")
plt.ylabel("Trace index (run order)")
plt.title(
    f"Peak Map (Background Subtracted)\n"
    f"Runs {valid_run_ids[0]} to {valid_run_ids[-1]} (DB order)"
)
plt.tight_layout()
plt.show()

# ---------------------
# SAVE DATA
# ---------------------
output_dir = "extracted_data"
os.makedirs(output_dir, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
base = os.path.join(output_dir, f"psd_map_{timestamp}")

# Binary (compact)
np.save(f"{base}_spectra.npy", spectra)
np.save(f"{base}_vgates.npy", vgates)
np.save(f"{base}_freqs.npy", freqs)
np.save(f"{base}_run_ids.npy", valid_run_ids)

# Text (human-readable)
np.savetxt(
    f"{base}_spectra.txt",
    spectra,
    delimiter="\t",
    header=(
        "PSD spectra matrix: rows = traces (DB order), cols = frequencies\n"
        f"Run IDs: {valid_run_ids[0]} to {valid_run_ids[-1]}"
    ),
)
np.savetxt(
    f"{base}_vgates.txt",
    vgates,
    delimiter="\t",
    header="Gate voltages V_gcs (mV) in same order as spectra rows",
)
np.savetxt(
    f"{base}_freqs.txt",
    freqs,
    delimiter="\t",
    header="Frequencies (MHz)",
)
np.savetxt(
    f"{base}_run_ids.txt",
    valid_run_ids,
    fmt="%d",
    header="Run IDs in same order as spectra rows",
)

with open(f"{base}_metadata.txt", "w") as f:
    f.write(f"Background ID: {background_id}\n")
    f.write(f"Run IDs requested: {run_ids[0]} to {run_ids[-1]}\n")
    f.write(f"Run IDs used (valid): {valid_run_ids[0]} to {valid_run_ids[-1]}\n")
    f.write(f"Number of spectra: {len(valid_run_ids)}\n")
    f.write(f"Number of frequency points: {len(freqs)}\n")
    f.write(f"Frequency range: {freqs[0]:.3f} - {freqs[-1]:.3f} MHz\n")
    f.write(f"Database: {qc.config['core']['db_location']}\n")

print(f"\n✅ Data saved in '{output_dir}' with timestamp {timestamp}")
print(f"   - {os.path.basename(base)}_spectra.npy/txt: {spectra.shape} matrix")
print(f"   - {os.path.basename(base)}_vgates.txt: {len(vgates)} values")
print(f"   - {os.path.basename(base)}_freqs.txt: {len(freqs)} values")
print(f"   - {os.path.basename(base)}_run_ids.txt: {len(valid_run_ids)} values")
print(f"   - {os.path.basename(base)}_metadata.txt: experiment info")
