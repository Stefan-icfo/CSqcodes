import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re
from matplotlib.colors import LogNorm

#qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v4.db'#for meas with other sideband
#qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v5.db'#for meas with lower sl drive
qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v7.db'#for meas on other side of cb peak and all measurements after diaster


#run_ids= list(range(29, 117, 2)) #for meas with other sideband,v4
#run_ids= list(range(66, 153, 2)) #for meas with lower SL drive,v5
#run_ids= list(range(96, 111, 2)) #before disaster, dropped due to hemt malfunction, v7; see nothing
#run_ids= list(range(161, 165, 2)) #before disaster, dropped due to hemt malfunction, v7; see nothing
#run_ids= list(range(337, 350, 2)) #after disaster, quick check if nonlinearities still exist, v7
#run_ids= list(range(357, 370, 2)) #after disaster, demod 10 Mhz detuned from mode, v7
#run_ids= list(range(376, 389, 2)) #after disaster, NOT DETUNED (mislabel), v7
run_ids= list(range(397, 410, 2)) #after disaster, detuned 10MHz, v7

data_name='avg_avg_psd_nodrive'
setpoint_name="freq_param"

freq,test_spectrum=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)
spectra = np.zeros((len(test_spectrum), len(run_ids)))

drives=[]
sum_spectra=[]

for i, run_id in enumerate(run_ids):
   exp_name,_, spectrum = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
   #print(f"extracted_runid {run_id}")
   spectra[:, i] = spectrum
   drive=get_metadata(run_id-1,print_it=False,return_data=True)['self_sigouts_sigouts1_amplitudes_amplitudes1_value']
   drives.append(drive*1e6)
   sum_spectra.append(spectrum.sum())


plt.plot(drives,sum_spectra,marker='o', linestyle='-',)
#plt.title('integral psd vs drive power; v5 rid 66-153')
#plt.title('integral psd vs drive power; v7 rid 337-349')
#plt.title('integral psd vs drive power; v7 rid 96-110')
#plt.title('integral psd vs drive power; v7 rid 161-165')
#plt.title('integral psd vs drive power; v7 rid 357-370')
#plt.title('integral psd vs drive power; v7 rid 376-388')
plt.title('integral psd vs drive power; v7 rid 397-409')
plt.xscale('log')
plt.xlabel('drive amp uVrms@instr')
plt.ylabel('integral psd [W]')
plt.show()
#plt.pcolormesh(freq[::-1] / 1e6, drives, spectra.T, cmap='magma', shading='auto')  # reverse frequency axis for other sideband

mesh = plt.pcolormesh(freq / 1e6, drives, spectra.T, cmap='magma', shading='auto',norm=LogNorm()) #log scale for post-disaster plots
#mesh = plt.pcolormesh(freq / 1e6, drives, spectra.T, cmap='magma', shading='auto')

# Set color limits on the returned QuadMesh object
#mesh.set_clim(0,1e-13)
# Use shading='auto' for better behavior
plt.yscale("log")  # Set Y-axis (drives) to log scale

plt.xlabel("Frequency (MHz)")
plt.ylabel("Drive Amplitude [uV]")
plt.colorbar(label="Amplitude [W/Hz]")
#plt.title("Log-scaled Drive vs Frequency Spectra;v7 rid 337-349")
#plt.title("Log-scaled Drive vs Frequency Spectra;v7 rid 96-110")
#plt.title("Log-scaled Drive vs Frequency Spectra;v7 rid 161-165")
#plt.title("Log-scaled Drive vs Frequency Spectra;v7 rid 357-370 - offres")
#plt.title("Log-scaled Drive vs Frequency Spectra;v7 rid 376-388 - onres")
plt.title("Log-scaled Drive vs Frequency Spectra;v7 rid 397-409 - offres")

plt.show()