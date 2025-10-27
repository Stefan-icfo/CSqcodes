import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import re

# ----------------------------------------
# CONFIG
# ----------------------------------------
qc.config["core"]["db_location"] = (
    r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v5.db"
)

qc.config["core"]["db_location"] = (
        r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v8.db")

background_id = 157
#all_data_ids = list(range(316, 358, 2))  # 316, 318, ..., 356
all_data_ids = list(range(66, 153, 2))  # 316, 318, ..., 356



background_id = 1158
#all_data_ids = list(range(316, 358, 2))  # 316, 318, ..., 356
all_data_ids = list(range(491, 561, 2))  # 316, 318, ..., 356


#background_id = 1158  # <<-- background run_id
#all_data_ids = list(range(1074, 1154, 2))  # g3 hole dot


#background_id = 1792#362
#all_data_ids = list(range(316, 358, 2))  # 316, 318, ..., 356
#background_id = 1469
#all_data_ids = list(range(568, 638, 2))  # 316, 318, ..., 356
#all_data_ids = list(range(1381, 1465, 2))  # 316, 318, ..., 356
#all_data_ids = list(range(1708, 1787, 2))  # 316, 318, ..., 356
#all_data_ids = list(range(1381, 1465, 2))  # 316, 318, ..., 356
#background_id = 1792
#all_data_ids = list(range(1708, 1787, 2))  # 316, 318, ..., 356

all_data_ids = list(range(1074, 1154, 2))  # 316, 318, ..., 356
#"""
#qc.config["core"]["db_location"] = (
#        r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v9.db")
#background_id = 314#1144
#all_data_ids = list(range(225, 309, 2))  # 316, 318, ..., 356
#all_data_ids = list(range(426, 508, 2))  # 316, 318, ..., 356
#all_data_ids = list(range(687, 750, 2))  # 316, 318, ..., 356
#all_data_ids = list(range(1077, 1140, 2))  # 316, 318, ..., 356
#"""
# ----------------------------------------
# ----------------------------------------
# HELPERS
# ----------------------------------------
def find_psd_key(param_data: dict) -> str:
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
# LOAD BACKGROUND
# ----------------------------------------
x_bg, y_bg, _ = load_psd(background_id)

# ----------------------------------------
# LOAD AND SUBTRACT
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

        # Just subtract point by point (truncate to min length)
        minlen = min(len(y), len(y_bg))
        y_corr = y[:minlen] - y_bg[:minlen]
        x_corr = x[:minlen]

        if x_common is None:
            x_common = x_corr

        # If lengths differ from x_common, truncate again
        minlen2 = min(len(x_common), len(y_corr))
        y_corr = y_corr[:minlen2]
        x_corr = x_corr[:minlen2]

        spectra.append(y_corr)
        drives.append(drive)
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


freq=x_common
freq_mask = (freq >= 151.175) & (freq <= 151.225) 


sum_spectra=np.sum(spectra,axis=1)
thermal_sum_mean=np.mean(sum_spectra[0:2])
print(f"mean thermal power {thermal_sum_mean}")
mean_spectra_reduced = np.sum(spectra[:,freq_mask], axis=1)

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
plt.clim(0,2e-6)
plt.xlabel("Frequency (MHz)", fontsize=12)
plt.ylabel("Drive amplitude (µV)", fontsize=12)
plt.title(f"Background-subtracted spectra runids {all_data_ids[0]} to {all_data_ids[-1]} in "+qc.config["core"]["db_location"][-15:-3], fontsize=14)
plt.colorbar(label="PSD (W/Hz)")
plt.tight_layout()
plt.show()




#"""
#extent = [freq[freq_mask][0], freq[freq_mask][-1], drives[0], drives[-1]]
#plt.imshow(
#    spectra[:,freq_mask],
#    aspect="auto",
#    extent=extent,
#    origin="lower",
#    cmap="plasma"
#)
#plt.clim(0,2e-6)
#plt.xlabel("Frequency (MHz)", fontsize=12)
#plt.ylabel("Drive amplitude (µV)", fontsize=12)
#plt.title("Background-subtracted spectra (point-by-point)", fontsize=14)
#plt.colorbar(label="PSD (W/Hz)")
#plt.tight_layout()
#plt.show()
#"""
print(f"✅ Loaded {len(valid_runs)} runs: {valid_runs}")
print(f"Drive amplitudes (µV): {drives}")

plt.plot(drives,sum_spectra/thermal_sum_mean,'*')
plt.plot(drives,mean_spectra_reduced/thermal_sum_mean,'go')
plt.xlabel("Drive amplitude (µV)", fontsize=12)
plt.ylabel("integrated spectra / montion amplitude squared", fontsize=12)
plt.xscale('log')
plt.show()

