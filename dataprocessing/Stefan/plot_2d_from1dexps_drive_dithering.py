import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re



data_name='avg_avg_psd_nodrive'
setpoint_name="freq_param"

#use background from dbv3
qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v3.db'

background_id=55
_,background_spectrum=extract_1d(background_id, data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)

print(f"background_nr_pt={len(background_spectrum)}")
qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v5.db'



#all_run_ids= list(range(274, 317, 2)) #for +mixdown
#excluded_ids ={290}# 

all_run_ids= list(range(318, 355, 2)) #for -mixdown
excluded_ids ={}#


run_ids = [rid for rid in all_run_ids if rid not in excluded_ids]



freq,test_spectrum=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)
spectra = np.zeros((len(test_spectrum), len(run_ids)))
print(f"data_nr_pt={len(test_spectrum)}")

drive_fs=[]

for i, run_id in enumerate(run_ids):
   exp_name,_, spectrum = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
   #print(f"extracted_runid {run_id}")
   spectra[:, i] = spectrum
   drive_f=get_metadata(run_id-1,print_it=False,return_data=True)['zurich_freq1']
   drive_fs.append(drive_f/1e6)

#plt.pcolormesh(freq[::-1] / 1e6, drives, spectra.T, cmap='magma', shading='auto')  # reverse frequency axis for other sideband
mesh = plt.pcolormesh(freq / 1e6,drive_fs,spectra.T,cmap='magma',shading='auto')

# Set color limits on the returned QuadMesh object
mesh.set_clim(0, 5e-14)

plt.xlabel("Frequency (MHz)")
plt.ylabel("Drive Frequency [MHz]")
plt.colorbar(label="Amplitude [W/Hz]")
plt.title("spectra vs drive_frequency,-sb,(v5;runs 317-356)")

plt.show()