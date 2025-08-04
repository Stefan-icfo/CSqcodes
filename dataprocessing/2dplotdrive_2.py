import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import re

# ------------------------------
# Database and run ID settings
# ------------------------------
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v2.db"
)
background_id = 1187
data_ids = list(range(951, 1216))
excluded_ids = [1184, 1185, 1186]
valid_ids = [i for i in data_ids if i not in excluded_ids]

# ------------------------------
# Helper functions
# ------------------------------
def load_psd(run_id):
    ds = qc.load_by_id(run_id)
    param_data = ds.get_parameter_data()
    key = list(param_data.keys())[0]
    y = param_data[key][key].flatten()
    x = param_data[key]["freq_param"].flatten() / 1e6
    sorted_idx = np.argsort(x)
    return x[sorted_idx], y[sorted_idx], ds

def extract_drive_amplitude_and_rep(exp_name):
    match = re.search(r"drive(\d+)u_rep(\d+)", exp_name)
    if match:
        return int(match.group(1)), int(match.group(2))
    return None, None

# ------------------------------
# Load background
# ------------------------------
x_bg, y_bg = load_psd(background_id)[:2]

# ------------------------------
# Load data
# ------------------------------
reps = {1: [], 2: []}
drive_amplitudes = {1: [], 2: []}
lengths = []

for run_id in valid_ids:
    try:
        x, y, ds = load_psd(run_id)
        drive, rep = extract_drive_amplitude_and_rep(ds.exp_name)
        if drive is None or rep not in (1, 2):
            continue
        minlen = min(len(y), len(y_bg))
        y_corr = y[:minlen] - y_bg[:minlen]
        reps[rep].append(y_corr)
        drive_amplitudes[rep].append(drive)
        lengths.append(minlen)
    except Exception as e:
        print(f"⚠️ Skipping {run_id}: {e}")

if not reps[1] or not reps[2]:
    print("❌ No valid spectra found for one or both repetitions.")
    print(f"Repetition 1 count: {len(reps[1])}")
    print(f"Repetition 2 count: {len(reps[2])}")
    raise RuntimeError("Aborting due to missing data.")

# ------------------------------
# Process and Plot
# ------------------------------
min_len = min(lengths)
x = x[:min_len]

rep1 = np.array([s[:min_len] for s in reps[1]])
rep2 = np.array([s[:min_len] for s in reps[2]])
mean_spectra = (rep1 + rep2) / 2

drive1 = np.array(drive_amplitudes[1])
drive2 = np.array(drive_amplitudes[2])
idx1 = np.argsort(drive1)
idx2 = np.argsort(drive2)
idxm = np.argsort((drive1 + drive2) / 2)

rep1 = rep1[idx1]
rep2 = rep2[idx2]
mean_spectra = mean_spectra[idxm]
drive1 = drive1[idx1]
drive2 = drive2[idx2]
drive_m = (drive1 + drive2) / 2

extent1 = [x[0], x[-1], drive1[0], drive1[-1]]
extent2 = [x[0], x[-1], drive2[0], drive2[-1]]
extentm = [x[0], x[-1], drive_m[0], drive_m[-1]]

# ------------------------------
# Plot all three maps
# ------------------------------
for data, extent, title in zip(
    [rep1, rep2, mean_spectra],
    [extent1, extent2, extentm],
    ["Repetition 1", "Repetition 2", "Mean of Reps"]
):
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        data,
        aspect="auto",
        extent=extent,
        origin="lower",
        cmap="plasma"
    )
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("Frequency (MHz)", fontsize=12)
    ax.set_ylabel("Drive Amplitude (μV)", fontsize=12)
    plt.colorbar(im, ax=ax, label="Amplitude (W/Hz)")
    plt.tight_layout()
    plt.show()

