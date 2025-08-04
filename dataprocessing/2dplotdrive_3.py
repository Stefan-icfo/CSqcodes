import numpy as np 
import matplotlib.pyplot as plt
import qcodes as qc
import re

# --- DATABASE CONFIGURATION ---
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v2.db"
)

background_id = 1187
excluded_ids = {1184, 1185, 1186, background_id, 1067, 1071,1075}

all_ids = list(range(1219, 1367, 2))
rep2_ids = [rid for rid in all_ids if rid not in excluded_ids]

# --- HELPER FUNCTIONS ---
def load_psd(run_id):
    ds = qc.load_by_id(run_id)
    pd = ds.get_parameter_data()
    key = list(pd.keys())[0]
    y = pd[key][key].flatten()
    x = pd[key]["freq_param"].flatten() / 1e6
    idx = np.argsort(x)
    return x[idx], y[idx], ds

def extract_drive(exp_name):
    m = re.search(r"drive(\d+)u", exp_name)
    return int(m.group(1)) if m else None

def is_rep2(name):
    return "rep2" in name.lower()

# --- LOAD BACKGROUND ---
x_bg, y_bg, _ = load_psd(background_id)

# --- LOAD REP2 DATA ONLY ---
spectra, drives, lengths = [], [], []
for run_id in rep2_ids:
    try:
        x, y, ds = load_psd(run_id)
        if not is_rep2(ds.exp_name):
            continue
        d = extract_drive(ds.exp_name)
        if d is None:
            continue
        L = min(len(y), len(y_bg))
        x_cut = x[:L]
        y_bg_interp = np.interp(x_cut, x_bg[:len(y_bg)], y_bg[:len(y_bg)])
        y_corr = y[:L] - y_bg_interp
        spectra.append(y_corr)
        drives.append(d)
        lengths.append(L)
        print(f"✔️ Loaded run {run_id} @ {d} μV")
    except Exception as e:
        print(f"⚠️ Skipping run {run_id}: {e}")

if not spectra:
    raise RuntimeError("❌ No valid rep2 data loaded.")

# --- FORMAT DATA FOR 2D GRID ---
min_len = min(lengths)
x = x_cut[:min_len]
spectra = np.array([s[:min_len] for s in spectra])
drives = np.array(drives)
sort_idx = np.argsort(drives)
spectra = spectra[sort_idx]
drives = drives[sort_idx]

# --- CREATE LOG SCALE Y GRID ---
log_drives = np.log10(drives)
log_drive_edges = np.concatenate([
    [log_drives[0] - (log_drives[1] - log_drives[0]) / 2],
    (log_drives[:-1] + log_drives[1:]) / 2,
    [log_drives[-1] + (log_drives[-1] - log_drives[-2]) / 2]
])
x_edges = np.concatenate([x, [x[-1] + (x[-1] - x[-2])]])
X, Y = np.meshgrid(x_edges, 10 ** log_drive_edges)

Z = spectra

# --- PLOT ---
plt.figure(figsize=(12, 6))
pc = plt.pcolormesh(X, Y, Z, cmap='magma', shading='auto')
plt.yscale('log')

plt.xlabel("Frequency (MHz)", fontsize=12)
plt.ylabel("Drive Amplitude (μV) [log scale]", fontsize=12)
plt.title("Rep 2 Spectra (Log Y Scale, Background Subtracted)", fontsize=14)
plt.colorbar(pc, label="PSD (W/Hz)")
plt.tight_layout()
plt.show()
