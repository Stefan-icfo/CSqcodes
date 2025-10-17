import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import re

# ----------------------------------------
# CONFIGURATION
# ----------------------------------------
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v14.db"
)

background_id = 295
#excluded_ids = [1184, 1185, 1186]
odd_data_ids = list(range(316, 358, 2))
all_data_ids =  odd_data_ids

# ----------------------------------------
# HELPERS
# ----------------------------------------
def load_psd(run_id):
    ds = qc.load_by_id(run_id)
    param_data = ds.get_parameter_data()
    key = list(param_data.keys())[0]
    y = param_data[key][key].flatten()
    x = param_data[key]["freq_param"].flatten() / 1e6
    sorted_idx = np.argsort(x)
    return x[sorted_idx], y[sorted_idx], ds

def extract_drive_amplitude(exp_name):
    match = re.search(r"drive(\d+)u", exp_name)
    if match:
        return int(match.group(1))  # in μV
    return None

# ----------------------------------------
# LOAD BACKGROUND
# ----------------------------------------
x_bg, y_bg = load_psd(background_id)[:2]

# ----------------------------------------
# LOAD DATA
# ----------------------------------------
spectra = []
drive_amplitudes = []
lengths = []

for run_id in all_data_ids:
    try:
        x, y, ds = load_psd(run_id)
        drive = extract_drive_amplitude(ds.exp_name)
        if drive is None:
            print(f"⚠️ Run {run_id}: Cannot extract drive from name.")
            continue
        minlen = min(len(y), len(y_bg))
        y_corr = y[:minlen] - y_bg[:minlen]
        spectra.append(y_corr)
        drive_amplitudes.append(drive)
        lengths.append(minlen)
    except Exception as e:
        print(f"⚠️ Skipping run {run_id}: {e}")
        continue

if len(spectra) == 0:
    raise RuntimeError("❌ No valid spectra loaded.")

# ----------------------------------------
# FORMAT DATA FOR IMSHOW
# ----------------------------------------
min_len = min(lengths)
spectra = np.array([s[:min_len] for s in spectra])
drive_amplitudes = np.array(drive_amplitudes)
x = x[:min_len]

# Sort by drive
sort_idx = np.argsort(drive_amplitudes)
spectra = spectra[sort_idx]
drive_amplitudes = drive_amplitudes[sort_idx]

# ----------------------------------------
# PLOT WITH IMSHOW
# ----------------------------------------
plt.figure(figsize=(12, 6))
extent = [x[0], x[-1], drive_amplitudes[0], drive_amplitudes[-1]]

plt.imshow(
    spectra,
    aspect="auto",
    extent=extent,
    origin="lower",
    cmap="plasma"
)

plt.xlabel("Frequency (MHz)", fontsize=12)
plt.ylabel("Drive Amplitude (μV)", fontsize=12)
plt.title("Background-subtracted Spectra", fontsize=14)
plt.colorbar(label="Amplitude (W/Hz)")
plt.tight_layout()
plt.show()
