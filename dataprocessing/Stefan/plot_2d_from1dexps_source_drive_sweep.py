import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re

qc.config["core"]["db_location"]=DATABASE_LOCATION

def extract_drive(exp_name):
    m = re.search(r"drive(\d+)u", exp_name)
    return int(m.group(1)) if m else None

run_ids= list(range(29, 117, 2)) 

data_name='avg_avg_psd_nodrive'
setpoint_name="freq_param"

freq,test_spectrum=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)
spectra = np.zeros((len(test_spectrum), len(run_ids)))

drives=[]

for i, run_id in enumerate(run_ids):
   exp_name,_, spectrum = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
   spectra[:, i] = spectrum
   drive=get_metadata(run_id-1,print=False,return_data=True)['self_sigouts_sigouts1_amplitudes_amplitudes1_value']
   drives.append(drive*1e6)

plt.pcolormesh(freq[::-1] / 1e6, drives, spectra.T, cmap='magma', shading='auto')  # Use shading='auto' for better behavior
plt.yscale("log")  # Set Y-axis (drives) to log scale

plt.xlabel("Frequency (MHz)")
plt.ylabel("Drive Amplitude [uV]")
plt.colorbar(label="Amplitude [W/Hz]")
plt.title("Log-scaled Drive vs Frequency Spectra")

plt.show()