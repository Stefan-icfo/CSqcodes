import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re
from matplotlib.colors import LogNorm

qc.config["core"]["db_location"] = r"C:\Users\sforstner\Desktop\Triton database\CD12_B5_F4v11.db"
print("Opening DB:", qc.config["core"]["db_location"])
#run_ids=[3263,3279,3295,3311,3327,3343,3359,3375,3391,3407,3423,3439,3455,3471,3487,3503]#
run_ids=[3554,3570,3602,3618,3263,3279,3295,3311,3327,3343,3359,3375,3391,3407,3423,3439,3455,3471,3487]#for best sensitivity ref linesweep 1946
run_ids=[3618,3263,3279,3295,3311,3327,3343,3359,3375,3391,3407,3423,3439,3455,3471,3487,3503]#for best sensitivity ref linesweep 1946
#electron_nrs=[3618,15,14,13,12,3327,10,9,8,3391,6,5,3439,3,3471,3487,3503]
#run_ids=[3690,3704,3718,3732,3746,3760,3774]#for constant sensitivity 5pA ref linesw 1946
#run_ids=[3116,3132,3148,3164,3180,3196,3212,3228,3244]#for ref linesweep 3104 best sens

############now some exclusion of fishy spectra######################
#run_ids=[3618,3263,3279,3295,3311,3343,3359,3375,3391,3407,3423,3455,3471,3487]#for best sensitivity ref linesweep 1946
run_ids=[3263,3279,3295,3311,3343,3359,3375,3391,3407,3423,3455]
#electron_nrs=[15,14,13,12,10,9,8,7,6,5,3]
electron_nrs=[15,14,13,12,10,9,8,7,6,5,3]
#electron_nrs=[15,14,13,12,10,9,8,6,5,3]
e_nr=True

areas=[]
frequencies=[]
g2_voltages=[]
 
for run_id in run_ids:
    metadata_temp=get_metadata(run_id-1,print_it=False,return_data=True)
    area=metadata_temp['integral_over_substracted_psd']
    g2_voltage=metadata_temp['qdac_ch02_dc_constant_V']
    frequency=metadata_temp['center_freq']
    print(area)
    areas.append(area)
    g2_voltages.append(g2_voltage)
    frequencies.append(frequency)
 

plt.plot(g2_voltages,areas,'g*')
plt.title("areas_vs_g2_voltages with best sensitivity ref linesweep 1946")
plt.xlabel("g2voltage")
plt.ylabel("PSD area")
for i, (x, y,run_id) in enumerate(zip(g2_voltages, areas,run_ids)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
plt.show()

if e_nr==True:
    plt.plot(electron_nrs,areas,'g*')
    plt.title("areas_vs_e_nr with best sensitivity ref linesweep 1946")
    plt.xlabel("nr electrons")
    plt.ylabel("PSD area")
    plt.xlim(0,16)
    plt.ylim(0,7e-14)
    for i, (x, y,run_id) in enumerate(zip(electron_nrs, areas,run_ids)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
plt.show()

plt.plot(g2_voltages,frequencies,'g*')
plt.title("frequencies_vs_e_nr with best sensitivity ref linesweep 1946")
plt.xlabel("g2voltage")
plt.ylabel("frequency")
for i, (x, y,run_id) in enumerate(zip(g2_voltages, frequencies,run_ids)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
plt.show()

if e_nr==True:
    plt.plot(electron_nrs,frequencies,'g*')
    plt.title("frequencies_vs_e_nr with best sensitivity ref linesweep 1946")
    plt.xlabel("nr electrons")
    plt.ylabel("frequency")
    for i, (x, y,run_id) in enumerate(zip(electron_nrs, frequencies,run_ids)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
plt.show()
 