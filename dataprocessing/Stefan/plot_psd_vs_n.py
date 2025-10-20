import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

#qc.config["core"]["db_location"] = r"C:\Users\sforstner\Desktop\Triton database\CD12_B5_F4v11.db"
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
#electron_nrs=[15,14,13,12,10,9,8,7,6,5,3]
#electron_nrs=[15,14,13,12,10,9,8,6,5,3]
e_nr=True

run_ids=[225,238,251,264,280,294,307,349,363,376,389,402,418,432,458,487,501]#in dbv1171025

electron_nrs=[4,5,6,7,8,9,10,12,13,14,15,16,17,18,20,22,23]

run_ids1=[225,238,251,264,280,294,307,349,376,389,402,418,432,458,487,501]#in dbv1171025

electron_nrs1=[4,5,6,7,8,9,10,12,14,15,16,17,18,20,22,23]

run_ids2=[527,540,553,566,579,595,609,647,661,673,682,691,700,709,721,747]#in dbv1171025

electron_nrs2=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]

run_ids=[527,540,553,566,579,595,609,647,661,673,682,691,700,709,721,747,225,238,251,264,280,294,307,349,376,389,402,418,432,458,487,501]#in dbv1171025
electron_nrs=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,4,5,6,7,8,9,10,12,14,15,16,17,18,20,22,23]


areas=[]
frequencies=[]
g2_voltages=[]
sensitivities=[]
 

areas1=[]
frequencies1=[]
g2_voltages1=[]
sensitivities1=[]

for run_id in run_ids1:#all
    metadata_temp=get_metadata(run_id-1,print_it=False,return_data=True)
    area=metadata_temp['integral_over_substracted_psd']
    g2_voltage=metadata_temp['qdac_ch02_dc_constant_V']
    frequency=metadata_temp['center_freq']
    sensitivity=metadata_temp['I_sens_sit']
    print(area)
    areas1.append(area)
    g2_voltages1.append(g2_voltage)
    frequencies1.append(frequency)
    sensitivities1.append(sensitivity)


g2_voltages2=[]
areas2=[]
sensitivities2=[]
frequencies2=[]
for run_id in run_ids2:
        metadata_temp=get_metadata(run_id-1,print_it=False,return_data=True)
        area=metadata_temp['integral_over_substracted_psd']
        g2_voltage=metadata_temp['qdac_ch02_dc_constant_V']
        frequency=metadata_temp['center_freq']
        sensitivity=metadata_temp['I_sens_sit']
        print(area)
        areas2.append(area)
        g2_voltages2.append(g2_voltage)
        frequencies2.append(frequency)
        sensitivities2.append(sensitivity)
    
 

plt.plot(g2_voltages1,areas1,'g*')
plt.plot(g2_voltages2,areas2,'r*')
plt.title("areas_vs_g2_voltages with best sensitivity ref linesweep 1946")
plt.xlabel("g2voltage")
plt.ylabel("PSD area")
for i, (x, y,run_id) in enumerate(zip(g2_voltages, areas,run_ids)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
plt.show()

if e_nr==True:
    plt.plot(electron_nrs1,areas1,'g*')
    plt.plot(electron_nrs2,areas2,'r*')
    plt.title("areas_vs_e_nr with best sensitivity ref linesweep 1946")
    plt.xlabel("nr electrons")
    plt.ylabel("PSD area")
    plt.xlim(left=0)
    plt.ylim(bottom=0) 
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for i, (x, y,run_id) in enumerate(zip(electron_nrs, areas,run_ids)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
    plt.show()

plt.plot(g2_voltages1,frequencies1,'g*')
plt.plot(g2_voltages2,frequencies2,'r*')
plt.title("frequencies_vs_e_nr with best sensitivity ref linesweep 1946")
plt.xlabel("g2voltage")
plt.ylabel("frequency")

for i, (x, y,run_id) in enumerate(zip(g2_voltages, frequencies,run_ids)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
plt.show()

if e_nr==True:
    plt.plot(electron_nrs1,frequencies1,'g*')
    plt.plot(electron_nrs2,frequencies2,'r*')
    plt.title("frequencies_vs_e_nr with best sensitivity ref linesweep 1946")
    plt.xlabel("nr electrons")
    plt.ylabel("frequency")
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for i, (x, y,run_id) in enumerate(zip(electron_nrs, frequencies,run_ids)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
    plt.show()
 
if e_nr==True:
    plt.plot(electron_nrs1,sensitivities1,'g*')
    plt.plot(electron_nrs2,sensitivities2,'r*')
    plt.title("sensitivities")
    plt.xlabel("nr electrons")
    plt.ylabel("sensitivities")
    plt.xlim(left=0)
    plt.ylim(bottom=0) 
    for i, (x, y,run_id) in enumerate(zip(electron_nrs, areas,run_ids)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
    plt.show()

    #plt.plot(electron_nrs,np.array(areas)*np.array(sensitivities)**2,'g*')
    #plt.title("scaled area")
    #plt.xlabel("nr electrons")
    #plt.ylabel("scaled psd")
    #plt.xlim(left=0)
    #plt.ylim(bottom=0) 
    #for i, (x, y,run_id) in enumerate(zip(electron_nrs, areas,run_ids)):
    #    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
    #plt.show()

    #plt.plot(electron_nrs,np.sqrt(areas)*np.array(sensitivities),'g*')
    #plt.title("scaled sqrt areas")
    #plt.xlabel("nr electrons")
    #plt.ylabel("scaled psd")
    #plt.xlim(left=0)
    #plt.ylim(bottom=0) 
    #for i, (x, y,run_id) in enumerate(zip(electron_nrs, areas,run_ids)):
    #    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
    #plt.show()

    plt.plot(electron_nrs1, np.sqrt(areas1)*np.array(sensitivities1), 'g*')
    plt.plot(electron_nrs2, np.sqrt(areas2)*np.array(sensitivities2), 'r*')

    plt.title("scaled sqrt areas")
    plt.xlabel("nr electrons")
    plt.ylabel("scaled psd")
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    # Calculate linear regression through origin first dataset
    x1 = np.array(electron_nrs1)
    y1 = np.sqrt(areas1) * np.array(sensitivities1)
    # For regression through origin: slope = sum(x*y) / sum(x^2)
    slope = np.sum(x1 * y1) / np.sum(x1**2)
    # Create regression line
    x1_line = np.linspace(0, max(electron_nrs), 100)
    y1_line = slope * x1_line
    # Plot the regression line
    plt.plot(x1_line, y1_line, 'g-', label=f'y = {slope:.4e}x', linewidth=2)
    # Add labels
    for i, (x, y, run_id) in enumerate(zip(electron_nrs1, np.sqrt(areas1)*np.array(sensitivities1), run_ids1)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
    
    # Calculate linear regression through origin first dataset
    x2 = np.array(electron_nrs2)
    y2 = np.sqrt(areas2) * np.array(sensitivities2)
    # For regression through origin: slope = sum(x*y) / sum(x^2)
    slope = np.sum(x2 * y2) / np.sum(x2**2)
    # Create regression line
    x2_line = np.linspace(0, max(electron_nrs), 100)
    y2_line = slope * x2_line
    # Plot the regression line
    plt.plot(x2_line, y2_line, 'r-', label=f'y = {slope:.4e}x', linewidth=2)
    # Add labels
    for i, (x, y, run_id) in enumerate(zip(electron_nrs2, np.sqrt(areas2)*np.array(sensitivities2), run_ids1)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')

    plt.legend()
    plt.show()


    

    plt.plot(electron_nrs2, np.sqrt(areas2)*np.array(sensitivities2), 'g*')

    plt.title("scaled sqrt areas")
    plt.xlabel("nr electrons")
    plt.ylabel("scaled psd")
    plt.xlim(left=0)
    plt.ylim(bottom=0)

  