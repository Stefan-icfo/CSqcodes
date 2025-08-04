import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import re

# --- CONFIGURATION ---
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v2.db"
)

background_id = 1187
run_ids = list(range(1375, 1423, 2))  # Inclusive of 1413

# --- HELPERS ---
def load_psd(run_id):
    ds = qc.load_by_id(run_id)
    param_data = ds.get_parameter_data()
    key = list(param_data.keys())[0]
    y = param_data[key][key].flatten()
    x = param_data[key]["freq_param"].flatten() / 1e6
    idx = np.argsort(x)
    return x[idx], y[idx], ds

def extract_drive(exp_name):
    m = re.search(r"_(\d+)u", exp_name)
    return int(m.group(1)) if m else None

# --- LOAD BACKGROUND ---
x_bg, y_bg, _ = load_psd(background_id)

# --- LOAD MEASUREMENTS ---
spectra = []
drives = []
lengths = []

for rid in run_ids:
    try:
        x, y, ds = load_psd(rid)
        drive = extract_drive(ds.exp_name)
        if drive is None:
            print(f"⚠️ Skipping {rid}: drive not found.")
            continue
        L = min(len(y), len(y_bg))
        x_cut = x[:L]
        y_bg_interp = np.interp(x_cut, x_bg[:len(y_bg)], y_bg[:len(y_bg)])
        y_corr = y[:L] - y_bg_interp
        spectra.append(y_corr)
        drives.append(drive)
        lengths.append(L)
        print(f"✔️ Loaded {rid} @ {drive} μV")
    except Exception as e:
        print(f"⚠️ Skipping {rid}: {e}")

if not spectra:
    raise RuntimeError("❌ No valid spectra loaded.")

# --- FORMAT DATA ---
min_len = min(lengths)
x = x_cut[:min_len]
spectra = np.array([s[:min_len] for s in spectra])
drives = np.array(drives)
sort_idx = np.argsort(drives)
spectra = spectra[sort_idx]
drives = drives[sort_idx]

# --- LOG SCALE Y GRID ---
log_drives = np.log10(drives)
log_drive_edges = np.concatenate([
    [log_drives[0] - (log_drives[1] - log_drives[0]) / 2],
    (log_drives[:-1] + log_drives[1:]) / 2,
    [log_drives[-1] + (log_drives[-1] - log_drives[-2]) / 2]
])
x_edges = np.concatenate([x, [x[-1] + (x[-1] - x[-2])]])
X, Y = np.meshgrid(x_edges, 10 ** log_drive_edges)

# --- PLOT ---
plt.figure(figsize=(12, 6))
pc = plt.pcolormesh(X, Y, spectra, cmap='magma', shading='auto')
plt.yscale('log')
plt.xlabel("Frequency (MHz)", fontsize=12)
plt.ylabel("Drive Amplitude (μV) [log scale]", fontsize=12)
plt.title("Spectra from run 1375 to 1413 (Background Subtracted)", fontsize=14)
plt.colorbar(pc, label="PSD (W/Hz)")
plt.tight_layout()
plt.show()


