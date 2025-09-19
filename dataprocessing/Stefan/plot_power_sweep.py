import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re







qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v8.db'#20mV run 10/09/25
qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v9.db'
data_name='I_rf'
setpoint_name="zurich_oscs0_freq"


run_ids= list(range(1299, 1307)) #v8 g2

run_ids= list(range(83, 94)) #v9 g4
run_ids= list(range(376, 417)) #v9 g4

#run_ids_cf= list(range(1162, 1169)) #v8 g3


linear_freq=141.8e6#g2
linear_freq=155.1935e6#g4
linear_freq_cf=151.2e6#g3

#freq,test_sweep=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)
#sweeps = np.zeros((len(test_sweep), len(run_ids)))

power_list=[]
maxf_list=[]



for i, run_id in enumerate(run_ids):
   exp_name,freq, sweep = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
   print(f"extracted_runid {run_id}")
   
   
   power=get_metadata(run_id,print_it=False,return_data=True)['self_sigouts_sigouts1_amplitudes_amplitudes1_value']
   power_list.append(power*1e3)
   plt.plot(freq,sweep,label=f'power={power*1e3:.3g} mV')
   Imax_id=np.argmax(sweep)
   f_max=freq[Imax_id]
   maxf_list.append((f_max-linear_freq)*1e3)

plt.xlabel("Frequency [MHz]")
plt.ylabel("drive@inst [mV]")
plt.legend()


plt.title("powersweep_g4dot")

plt.show()
"""
power_list_cf=[]
maxf_list_cf=[]

for i, run_id in enumerate(run_ids_cf):
   exp_name,freq, sweep = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
   print(f"extracted_runid {run_id}")
   
   
   power=get_metadata(run_id,print_it=False,return_data=True)['self_sigouts_sigouts1_amplitudes_amplitudes1_value']
   power_list_cf.append(power*1e3)
   plt.plot(freq,sweep,label=f'power={power*1e3:.3g} mV')
   Imax_id=np.argmax(sweep)
   f_max=freq[Imax_id]
   maxf_list_cf.append((f_max-linear_freq_cf)*1e3)
#plt.pcolormesh(freq[::-1] / 1e6, drives, spectra.T, cmap='magma', shading='auto')  # reverse frequency axis for other sideband




plt.xlabel("Frequency [kHz]")
plt.ylabel("drive@inst [mV]")
plt.legend()


plt.title("powersweep_g3dot")

plt.show()
"""



plt.plot(power_list,maxf_list,'*',label="g4dot")
#plt.plot(power_list_cf,maxf_list_cf,'*',label="g3dot")
plt.legend()
plt.xscale('log') 
plt.ylabel("Frequency [MHz]")
plt.xlabel("drive")
plt.title("duffing frequency drag vs power")
plt.show()
