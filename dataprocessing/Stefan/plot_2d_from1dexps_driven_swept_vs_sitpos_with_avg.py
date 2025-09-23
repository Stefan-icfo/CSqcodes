import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re







data_name='I_rf_avg'
setpoint_name="zurich_oscs0_freq"



qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v11.db'#20mV run 10/09/25

avg_num=10
excluded_ids ={}

all_run_ids= list(range(74, 124 ,avg_num)) #v9 100uV on gate2 at instr; g3 dot;runs from 11 &13 /09/25
excluded_ids ={}

run_ids = [rid for rid in all_run_ids if rid not in excluded_ids]



freq,test_sweep=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)
sweeps = np.zeros((len(test_sweep), len(run_ids)))

sitpos_list=[]

for i, run_id in enumerate(run_ids):
   #fill initial sweep
   exp_name,_, sweep = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
   print(f"extracted_runid {run_id} of first run at this Voltage")
   sweeps[:, i] = sweep
   for n in range(avg_num-1):
        exp_name,_, sweep = extract_1d(run_id+n+1, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
        print(f"extracted_runid {run_id}")
        sweeps[:, i] = sweeps[:, i]+sweep
   sweeps[:, i]=sweeps[:, i]/avg_num #normalize...remove if doesnt work
   sitpos=get_metadata(run_id-1,print_it=False,return_data=True)['qdac_ch06_dc_constant_V']
   sitpos_list.append(sitpos*1e3)

#plt.pcolormesh(freq[::-1] / 1e6, drives, spectra.T, cmap='magma', shading='auto')  # reverse frequency axis for other sideband
mesh = plt.pcolormesh(freq / 1e6,sitpos_list,sweeps.T*1e12,cmap='magma',shading='auto')

# Set color limits on the returned QuadMesh object
#mesh.set_clim(0,2)

plt.xlabel("Frequency [MHz]")
plt.ylabel("sitpos [mV]")
plt.colorbar(label="current [pA]")
#plt.title("sweeps vs sitpos,-sb,100uV@instrg2(v7;runs 14-44)")
#plt.title("sweeps vs sitpos,-sb,500uV@instrg2(v7;runs 244-284)")
#plt.title("sweeps vs sitpos,-sb,100uV@instrg2(v6;runs 26-43)")
#plt.title("sweeps vs sitpos,-sb,1mV@instrg2(v7;runs 285-325)")
#plt.title("sweeps vs sitpos,-sb,10mV@instrg2(v8;runs 141-171)")
#plt.title("sweeps vs sitpos,-sb,10mV@instrg2(v8;runs 690-730)")
#plt.title("sweeps vs sitpos,-sb,1.25mV@instrg2(v8;runs 772-812)")
#plt.title("sweeps vs sitpos,-sb,327uV@instrg2(v8;runs 953-993)")
plt.title("sweeps vs sitpos,-sb,312uV@instrg2(v9;runs 1261-1292)")
plt.title("sweeps vs sitpos,-sb,20mV@instrg2g4dot(v9;runs 97-137)")
plt.title("sweeps vs sitpos,-sb,5mV@instrg2g4dot(v9;runs 138-178)")
plt.title("sweeps vs sitpos,-sb,1.25mV@instrg2g4dot(v9;runs 179-219)")
plt.show()
