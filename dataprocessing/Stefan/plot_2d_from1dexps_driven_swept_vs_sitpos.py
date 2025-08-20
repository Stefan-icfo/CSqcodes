import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re





qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v7.db'#500uV run after disaster and 100uV run before disaster; 1mV run after disaster #runs on first cb peak after before collecting mixture 

#qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v6.db'#100uV run before disaster


data_name='I_rf'
setpoint_name="zurich_oscs0_freq"

#all_run_ids= list(range(26, 43)) #v6 100uVrf on gate @instr
#excluded_ids ={}#
 

#all_run_ids= list(range(285, 326)) #v7 1mV at gate on instr
#excluded_ids ={}#

#all_run_ids= list(range(14, 44)) #v7; for -mixdown; quick run before the disaster; 100uVrf on gate @instr
#excluded_ids ={}#

#all_run_ids= list(range(244, 284)) #v6; for -mixdown; after disaster (13/08/25); 500Vrf on gate@instr. using second peak
#excluded_ids ={}#

all_run_ids= list(range(515, 555)) #v7 2mV on gate at instr; before collecting mixture
excluded_ids ={}#

run_ids = [rid for rid in all_run_ids if rid not in excluded_ids]



freq,test_sweep=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)
sweeps = np.zeros((len(test_sweep), len(run_ids)))

sitpos_list=[]

for i, run_id in enumerate(run_ids):
   exp_name,_, sweep = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
   print(f"extracted_runid {run_id}")
   sweeps[:, i] = sweep
   sitpos=get_metadata(run_id-1,print_it=False,return_data=True)['qdac_ch06_dc_constant_V']
   sitpos_list.append(sitpos*1e3)

#plt.pcolormesh(freq[::-1] / 1e6, drives, spectra.T, cmap='magma', shading='auto')  # reverse frequency axis for other sideband
mesh = plt.pcolormesh(freq / 1e6,sitpos_list,sweeps.T*1e12,cmap='magma',shading='auto')

# Set color limits on the returned QuadMesh object
mesh.set_clim(0,5)

plt.xlabel("Frequency [MHz]")
plt.ylabel("sitpos [mV]")
plt.colorbar(label="current [pA]")
#plt.title("sweeps vs sitpos,-sb,100uV@instrg2(v7;runs 14-44)")
#plt.title("sweeps vs sitpos,-sb,500uV@instrg2(v7;runs 244-284)")
#plt.title("sweeps vs sitpos,-sb,100uV@instrg2(v6;runs 26-43)")
#plt.title("sweeps vs sitpos,-sb,1mV@instrg2(v7;runs 285-325)")
plt.title("sweeps vs sitpos,-sb,2mV@instrg2(v7;runs 514-556)")

plt.show()
