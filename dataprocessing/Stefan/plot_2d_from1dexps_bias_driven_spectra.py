import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re
from matplotlib.colors import LogNorm

qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v6.db'#testrun - seems only feasible one



run_ids= list(range(49, 66, 3)) #before disaster, testrun

data_name='avg_avg_psd_nodrive'
setpoint_name="freq_param"

freq,test_spectrum=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)
spectra = np.zeros((len(test_spectrum), len(run_ids)))

biases=[]
sum_spectra=[]

for i, run_id in enumerate(run_ids):
   exp_name,_, spectrum = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=True,return_exp_name=True)
   print(f"extracted_runid {run_id}")
   spectra[:, i] = spectrum
   bias=get_metadata(run_id-1,print_it=False,return_data=True)['qdac_ch07_dc_constant_V']
   biases.append(bias*1e6)
   sum_spectra.append(spectrum.sum())




mesh = plt.pcolormesh(freq / 1e6, biases, spectra.T, cmap='magma', shading='auto',norm=LogNorm()) #log scale for post-disaster plots
#mesh = plt.pcolormesh(freq / 1e6, drives, spectra.T, cmap='magma', shading='auto')

# Set color limits on the returned QuadMesh object
#mesh.set_clim(0,1e-13)
# Use shading='auto' for better behavior
#plt.yscale("log")  # Set Y-axis (drives) to log scale

plt.xlabel("Frequency (MHz)")
plt.ylabel("Drive Amplitude [uV]")
plt.colorbar(label="Amplitude [W/Hz]")
#plt.title("Log-scaled Drive vs Frequency Spectra;v7 rid 337-349")
#plt.title("Log-scaled Drive vs Frequency Spectra;v7 rid 96-110")
#plt.title("Log-scaled Drive vs Frequency Spectra;v7 rid 161-165")
#plt.title("Log-scaled Drive vs Frequency Spectra;v7 rid 357-370 - offres")
#plt.title("Log-scaled Drive vs Frequency Spectra;v7 rid 376-388 - onres")
plt.title("bias vs Frequency Spectra;v6 rid 49-64 (testrun but only working one)")

plt.show()