import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import re

# ---------------------
# CONFIGURATION
# ---------------------
#qc.config["core"]["db_location"] = (
    #"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v8.db"
#)
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v9.db"
)

#background_id = 1343
#run_ids = list(range(1345, 1473,2))

#background_id = 1343
#run_ids = list(range(1345, 1470,2))


#run_ids = list(range(1483,1581,2))
#background_id = 1585
#run_ids = list(range(1587, 1681,2))

#background_id = 1481
#run_ids = list(range(1483, 1582,2))

background_id = 1684
run_ids = list(range(1686, 1931,2))


#background_id = 790
#run_ids = list(range(792, 965,2))


background_id = 2610
run_ids = list(range(2612, 2868,2))


background_id = 3664
run_ids = list(range(3666, 3780,2))


qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v21_27_10_25.db"
)

background_id = 540
run_ids = list(range(542, 639,2))

run_ids = list(range(400, 524,2))
background_id = 252
run_ids = list(range(254, 378,2))
run_ids = list(range(108, 232,2))



#qc.config["core"]["db_location"] = (
#    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v24_01_11_25.db"
#)

background_id = 182
run_ids = list(range(184, 245,2))#6 electrons


#background_id = 415
#run_ids = list(range(417, 478,2))#11 electrons


#background_id = 648
#run_ids = list(range(650, 711,2))#16 electrons

#background_id = 881
#run_ids = list(range(883, 944,2))#21 electrons



#qc.config["core"]["db_location"] = (
#    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v25_02_11_25.db"
#)

#background_id = 148
#run_ids = list(range(150, 211,2))#26 electrons

#now more precise softening done on Saturday 151125
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v30_14_11_25.db"
)

background_id = 268
run_ids = list(range(270, 487,2))#I think 10 e, check!

background_id = 578
run_ids = list(range(580, 743,2))#I think 15 e, check!


qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v31_17_11_25.db"
)

background_id = 34
run_ids = list(range(36, 251,2))#I think 15 e, check!

#background_id = 768
#run_ids = list(range(770, 853,2))#I think 15 e, check!


background_id = 370
run_ids = list(range(372, 544,2))#I think 15 e, check!
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
#def extract_gcs_voltage(name):
#    match = re.search(r"quick=([\d.]+)", name)
#    return float(match.group(1)) if match else None
def extract_gcs_voltage(name):
    match = re.search(r"gcs_=([\d.]+)\s*mV", name)
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
print(lengths)
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
plt.imshow(spectra, aspect='auto', extent=extent, origin='lower', 
           cmap='viridis', interpolation='none')
plt.colorbar(label='Peak PSD Intensity (W/Hz)')
plt.xlabel("Peak Frequency (MHz)")
plt.ylabel("Gate Voltage V_gcs (mV)")
plt.title(f" Peak Map (Background Subtracted) runids {run_ids[0]} to {run_ids[-1]}")
plt.tight_layout()
plt.show()


# ---------------------
# SAVE DATA TO FILES
# ---------------------
import os
from datetime import datetime

# Create output directory if it doesn't exist
output_dir = "extracted_data"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Generate timestamp for unique filenames
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# Save as NumPy arrays (binary format - fast and compact)
np.save(f"{output_dir}/spectra_{timestamp}.npy", spectra)
np.save(f"{output_dir}/vgates_{timestamp}.npy", vgates)
np.save(f"{output_dir}/freqs_{timestamp}.npy", freqs)

# Also save as text files (human-readable)
np.savetxt(f"{output_dir}/spectra_{timestamp}.txt", spectra, delimiter='\t', 
           header=f"PSD spectra matrix: rows=gate voltages, cols=frequencies\nRun IDs: {run_ids[0]} to {run_ids[-1]}")
np.savetxt(f"{output_dir}/vgates_{timestamp}.txt", vgates, delimiter='\t',
           header="Gate voltages V_gcs (mV)")
np.savetxt(f"{output_dir}/freqs_{timestamp}.txt", freqs, delimiter='\t',
           header="Frequencies (MHz)")

# Save metadata
with open(f"{output_dir}/metadata_{timestamp}.txt", 'w') as f:
    f.write(f"Background ID: {background_id}\n")
    f.write(f"Run IDs: {run_ids[0]} to {run_ids[-1]}\n")
    f.write(f"Number of spectra: {len(vgates)}\n")
    f.write(f"Number of frequency points: {len(freqs)}\n")
    f.write(f"Frequency range: {freqs[0]:.3f} - {freqs[-1]:.3f} MHz\n")
    f.write(f"Gate voltage range: {vgates[0]:.3f} - {vgates[-1]:.3f} mV\n")
    f.write(f"Database: {qc.config['core']['db_location']}\n")

print(f"\n✅ Data saved to '{output_dir}/' with timestamp {timestamp}")
print(f"   - spectra_{timestamp}.npy/txt: {spectra.shape} matrix")
print(f"   - vgates_{timestamp}.txt: {len(vgates)} values")
print(f"   - freqs_{timestamp}.txt: {len(freqs)} values")
print(f"   - metadata_{timestamp}.txt: experiment info")




