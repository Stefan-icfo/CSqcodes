import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import re

# ----------------------------------------
# CONFIG
# ----------------------------------------
qc.config["core"]["db_location"] = (
    r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v8.db"
)

all_data_ids = list(range(316, 358, 2))  # 316, 318, ..., 356
all_data_ids = list(range(568, 638, 2))  # 316, 318, ..., 356
all_data_ids = list(range(491, 561, 2))  # 316, 318, ..., 356
all_data_ids = list(range(1074, 1154, 2))  # 316, 318, ..., 356
all_data_ids = list(range(1381, 1465, 2))  # 316, 318, ..., 356
all_data_ids = list(range(1708, 1787, 2))  # 316, 318, ..., 356

qc.config["core"]["db_location"] = (
    r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v9.db"
)

all_data_ids = list(range(225, 309, 2))  # 316, 318, ..., 356
# ----------------------------------------
# HELPERS
# ----------------------------------------
def find_psd_key(param_data: dict) -> str:
    """Find a key that contains 'freq_param'. Fallback: first key."""
    for k, v in param_data.items():
        if "freq_param" in v:
            return k
    return list(param_data.keys())[0]

def load_psd(run_id):
    ds = qc.load_by_id(run_id)
    p = ds.get_parameter_data()
    key = find_psd_key(p)
    y = np.asarray(p[key][key]).flatten()
    x = np.asarray(p[key]["freq_param"]).flatten() / 1e6  # MHz
    idx = np.argsort(x)
    return x[idx], y[idx], ds

def extract_drive_amplitude(exp_name: str):
    """Extract drive amplitude in µV from exp_name."""
    patterns = [
        r"[dD]rive\s*([0-9]+)\s*u[V]?",
        r"_(\d+)u_(?=1Dcs)",
        r"_(\d+)u(?:_|$)",
        r"\b(\d+)u\b",
    ]
    for pat in patterns:
        m = re.search(pat, exp_name)
        if m:
            return int(m.group(1))
    return None

# ----------------------------------------
# LOAD ALL DATA
# ----------------------------------------
spectra = []
drives = []
valid_runs = []
x_common = None

for run_id in all_data_ids:
    try:
        x, y, ds = load_psd(run_id)
        drive = extract_drive_amplitude(ds.exp_name)
        if drive is None:
            print(f"⚠️ Run {run_id}: could not extract drive from '{ds.exp_name}'")
            continue

        if x_common is None:
            x_common = x  # use the first run’s grid as reference

        # Interpolate to common grid if needed
        if len(x) != len(x_common) or not np.allclose(x, x_common):
            y_interp = np.interp(x_common, x, y, left=np.nan, right=np.nan)
        else:
            y_interp = y

        spectra.append(y_interp)
        drives.append(drive/1e3)
        valid_runs.append(run_id)

    except Exception as e:
        print(f"⚠️ Run {run_id}: skipped → {e}")

if len(spectra) == 0:
    raise RuntimeError("❌ No valid spectra loaded.")

spectra = np.array(spectra)
drives = np.array(drives)

# Sort by drive
order = np.argsort(drives)
spectra = spectra[order]
drives = drives[order]

# ----------------------------------------
# PLOT
# ----------------------------------------
plt.figure(figsize=(12, 6))
extent = [x_common[0], x_common[-1], drives[0], drives[-1]]

plt.imshow(
    spectra,
    aspect="auto",
    extent=extent,
    origin="lower",
    cmap="plasma"
)

plt.xlabel("Frequency (MHz)", fontsize=12)
plt.ylabel("Drive amplitude (mV)", fontsize=12)
plt.title("Raw spectra (no background subtraction)", fontsize=14)
plt.colorbar(label="PSD (a.u).")
plt.tight_layout()
plt.show()

#first_half = spectra[:len(spectra)//3]
avg_spectra = np.mean(spectra, axis=0)
plt.plot(x,avg_spectra)
plt.xlabel("Frequency (MHz)", fontsize=12)
plt.ylabel("psd", fontsize=12)
plt.title("Raw spectra (no background subtraction)", fontsize=14)

plt.show()

print(f"✅ Loaded {len(valid_runs)} runs: {valid_runs}")
print(f"Drive amplitudes (µV): {drives}")
