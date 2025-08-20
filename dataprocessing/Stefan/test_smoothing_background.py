import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re


qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v3.db'

def voltage_to_psd(v_rms, rbw, impedance=50):
  
    # Calculate PSD using the formula
    psd = (v_rms ** 2) / (impedance * rbw)
    return psd

#run_ids= list(range(29, 117, 2)) #for meas with other sideband
run_ids= list(range(55, 94, 2)) #for meas with lower SL drive# measurement of 4.8.
background_id=95

#data_name='avg_avg_psd_nodrive'
data_name='V_fft_avg_avg'
setpoint_name="freq_param"

freq,test_spectrum=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)

freq,V_background=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)
spectra = np.zeros((len(test_spectrum), len(run_ids)))
rbw=get_metadata(run_ids[0]-1,print_it=False,return_data=True)['rbw']#take the rbw of the first driven run for all conversions

probe_amps=[]

for i, run_id in enumerate(run_ids):
   _, V_spectrum = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False)
   V_spectrum=V_spectrum-V_background
   p_spectrum=voltage_to_psd(v_rms=V_spectrum,rbw=rbw)
   #print(f"extracted_runid {run_id}")
   #spectra[:, i] = spectrum
   probe_amp=get_metadata(run_id-1,print_it=False,return_data=True)['self_sigouts_sigouts0_amplitudes_amplitudes0_value']
   probe_amps.append(probe_amp*1e3)
   spectra[:, i] = p_spectrum/probe_amp**4#scale it py the probe amp squared

#plt.pcolormesh(freq[::-1] / 1e6, drives, spectra.T, cmap='magma', shading='auto')  # reverse frequency axis for other sideband
plt.pcolormesh(freq / 1e6, probe_amps, spectra.T, cmap='magma', shading='auto')  # Use shading='auto' for better behavior
#plt.yscale("log")  # Set Y-axis (drives) to log scale

plt.xlabel("Frequency (MHz)")
#plt.ylabel("Probe Amplitude [mV]")
plt.ylabel("Scaled Probe Amplitude [mV]")
plt.colorbar(label="Amplitude [W/Hz]")
plt.title("probe amp spectrum")

plt.show()