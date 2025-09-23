import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re








qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v9.db'
data_name='I_rf'
setpoint_name="zurich_oscs0_freq"


#run_ids= list(range(1299, 1307)) #v8 g2

run_ids= list(range(2028,2078)) #v9 g4


#run_ids_cf= list(range(1162, 1169)) #v8 g3


sweeps=[]

for i, run_id in enumerate(run_ids):
   exp_name,freq, sweep = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
   print(f"extracted_runid {run_id}")
   
   sweeps.append(sweep)
   
avg_sweep=np.mean(sweeps,axis=0)  
plt.plot(freq,avg_sweep)
plt.xlabel("Frequency [MHz]")
plt.ylabel("current[pA]")
#plt.legend()


plt.title("avg100uVdrive")

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



