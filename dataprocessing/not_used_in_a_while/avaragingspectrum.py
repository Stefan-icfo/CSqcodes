import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ----------- Load database and dataset -------------------------------
qc.config["core"]["db_location"] = "C:/Users/LAB-nanooptomechanic/Documents/MartaStefan/CSqcodes/Data/Raw_data/CD12_B5_F4v2.db"
dataset = qc.load_by_id(114)

# ----------- Extract data --------------------------------------------
data = dataset.get_parameter_data('avg_psd')['avg_psd']
freq = np.array(data['freq_param'])
time = np.array(data['time_param'])
psd_raw = np.array(data['avg_psd'])

# ----------- Reshape into 2D matrix (time x freq) --------------------
unique_time = np.unique(time)
unique_freq = np.unique(freq)
n_time = len(unique_time)
n_freq = len(unique_freq)

psd = psd_raw.reshape(n_time, n_freq)
freq = unique_freq
time = unique_time

# ----------- Mask based on time range to remove ----------------------
bad_mask = (time >= 4806) & (time <= 8150)
good_mask = ~bad_mask

# ----------- PLOT 1: All traces (black) and excluded (red) ----------
plt.figure(figsize=(10, 6))
for i in range(n_time):
    color = 'red' if bad_mask[i] else 'black'
    plt.plot(freq, psd[i], color=color, alpha=0.07, linewidth=0.5)
plt.xlabel("Frequency (Hz)")
plt.ylabel("PSD (W/Hz)")
plt.title("All traces (black) and excluded (red)")
plt.grid(True)
plt.tight_layout()
plt.show()

# ----------- PLOT 2: Mean with and without excluded -------------------
mean_all = np.mean(psd, axis=0)
mean_good = np.mean(psd[good_mask], axis=0)

plt.figure(figsize=(10, 6))
plt.plot(freq, mean_good, color='black', label="Mean (good only)")
plt.plot(freq, mean_all, color='red', linestyle='--', label="Mean (all)")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Averaged PSD (W/Hz)")
plt.title("Mean PSD: good vs all traces")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# ----------- PLOT 3: Only mean of good traces ------------------------
plt.figure(figsize=(10, 6))
plt.plot(freq, mean_good, color='black')
plt.xlabel("Frequency (Hz)")
plt.ylabel("Averaged PSD (W/Hz)")
plt.title("Mean PSD (only good traces)")
plt.grid(True)
plt.tight_layout()
plt.show()

# ----------- CENTERING AND FINAL AVERAGE -----------------------------
ref_y_index = 1000
psd_good = psd[good_mask]
ref_spectrum = psd_good[ref_y_index]
peak_index = np.argmax(ref_spectrum)

pad_width = max(peak_index, n_freq - peak_index - 1)
centered_psd = []

for row in psd_good:
    centered = np.zeros(2 * pad_width + 1)
    start = pad_width - peak_index
    centered[start:start + n_freq] = row
    centered_psd.append(centered)

centered_psd = np.array(centered_psd)
mean_centered_psd = np.mean(centered_psd, axis=0)

df = freq[1] - freq[0]
freq_centered = (np.arange(len(mean_centered_psd)) - pad_width) * df

# ----------- PLOT 4: Centered average ---------------------------------
plt.figure(figsize=(10, 6))
plt.plot(freq_centered + freq[peak_index], mean_centered_psd)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Averaged PSD (W/Hz)")
plt.title("Centered and averaged PSD (only good traces)")
plt.grid(True)
plt.tight_layout()
plt.show()
