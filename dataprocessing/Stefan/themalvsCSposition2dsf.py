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
#qc.config["core"]["db_location"] = (
#    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v30_14_11_25.db"
#)

#background_id = 268
#run_ids = list(range(270, 487,2))#I think 10 e, check!

#background_id = 578
#run_ids = list(range(580, 743,2))#I think 15 e, check!


#qc.config["core"]["db_location"] = (
#    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v31_17_11_25.db"
#)

#background_id = 34
#run_ids = list(range(36, 251,2))#I think 15 e, check!

#background_id = 768
#run_ids = list(range(770, 853,2))#I think 15 e, check!


#background_id = 370
#run_ids = list(range(372, 544,2))#I think 15 e, check!
#signal_key = "avg_avg_psd_nodrive"
#freq_key = "freq_param"



#qc.config["core"]["db_location"] = (
#    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v32_19_11_25.db"#20 and 21 11 - sweeping g1 in single e thermal 150M config
#)
#background_id = 500
#run_ids = list(range(503, 561,3))#
#signal_key = "avg_avg_psd_nodrive"
#freq_key = "freq_param"


#qc.config["core"]["db_location"] = (
#    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v31_17_11_25.db"#20 and 21 11 - sweeping g2 in single e thermal 150M config
#)
#background_id = 667
#run_ids = list(range(670, 768,3))#
#signal_key = "avg_avg_psd_nodrive"
#freq_key = "freq_param"


#qc.config["core"]["db_location"] = (
#    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v33_21_11_25.db"#26 and 21e, 221125,vgcs
#)
#background_id = 43
#run_ids = list(range(45, 118,2))#26e

#background_id = 128
#run_ids = list(range(130, 203,2))#21e

qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v33_21_11_25.db"#26 and 21e, 221125,vgcs
)

background_id = 529
run_ids = list(range(532, 587,3))#

background_id = 470
run_ids = list(range(473, 519,3))#


qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v38_28_11_25.db"#holes
)

background_id = 10
run_ids = list(range(12, 193,5))#1 hole

background_id = 203
run_ids = list(range(205, 350,4))#2 holes


background_id = 359
run_ids = list(range(361, 506,4))#6 holes

qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v39_29_11_25.db"#holes
)


background_id = 15
run_ids = list(range(17, 138,4))#10 holes

#background_id = 173
#run_ids = list(range(175, 321,4))#16 holes

#background_id = 331
#run_ids = list(range(333, 478,4))#20 holes


background_id = 489
run_ids = list(range(491, 636,4))#24 holes

#background_id = 645
#run_ids = list(range(647, 744,4))#


#qc.config["core"]["db_location"] = (
#    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v40_30_11_25.db"#holes
#)




#background_id = 7
#run_ids = list(range(9, 66,4))#12 holes

#background_id = 83
#run_ids = list(range(85, 230,4))#28 holes

#background_id = 240
#run_ids = list(range(242, 288,5))#6 holes, with readjustment


qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v41_01_12_25.db"#holes
)


background_id = 125
run_ids = list(range(127, 272,4))#30 holes, with readjustment

background_id = 284
run_ids = list(range(286, 411,4))#15 holes, with readjustment


qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v42_02_12_25.db"#holes
)

background_id = 184
run_ids = list(range(186, 331,4))#20 holes, with readjustment

background_id = 26
run_ids = list(range(28, 173,4))#20 holes, with readjustment


qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v46_07_12_25.db"#holes
)


background_id = 65
run_ids = list(range(67, 217,3))

background_id = 226
run_ids = list(range(228, 421,3))

#background_id = 425
#run_ids = list(range(427, 527,3))

background_id = 538
run_ids = list(range(540, 667,3))


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
plt.xlabel("Peak Frequency (MHz) ")
plt.ylabel("Gate Voltage V_gcs (mV)")
plt.title(f"Peak Map (Background Sub) ids {run_ids[0]} to {run_ids[-1]} in "+qc.config["core"]["db_location"][-15:-3])
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




