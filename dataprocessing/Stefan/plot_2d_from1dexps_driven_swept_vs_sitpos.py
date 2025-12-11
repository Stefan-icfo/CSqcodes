import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re





qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v7.db'#500uV run after disaster and 100uV run before disaster; 1mV run after disaster #runs on first cb peak after before collecting mixture 

#qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v6.db'#100uV run before disaster

qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v8.db'#20mV run 10/09/25

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

#all_run_ids= list(range(515, 555)) #v7 2mV on gate at instr; before collecting mixture
#excluded_ids ={}#

#all_run_ids= list(range(141, 171)) #v8 20mV on gate at instr; runs from 10/09/25
#excluded_ids ={}#

#all_run_ids= list(range(649, 689)) #v8 10mV on gate at instr; g3 dot;runs from 11 &13 /09/25
#excluded_ids ={}#

#all_run_ids= list(range(690, 730)) #v8 5mV on gate at instr; g3 dot;runs from 11 &13 /09/25
#excluded_ids ={}#

all_run_ids= list(range(731, 771)) #v8 2.5mV on gate at instr; g3 dot;runs from 11 &13 /09/25
excluded_ids ={}#

all_run_ids= list(range(772, 812)) #v8 1.25mV on gate at instr; g3 dot;runs from 11 &13 /09/25
excluded_ids ={}#

all_run_ids= list(range(856, 896)) #v8 625uV on gate at instr; g3 dot;runs from 11 &13 /09/25
excluded_ids ={}#




all_run_ids= list(range(953, 993)) #v8 625uV on gate at instr; g3 dot;runs from 11 &13 /09/25
excluded_ids ={}#

all_run_ids= list(range(1245, 1259)) #v8 625uV on gate at instr; g3 dot;runs from 11 &13 /09/25
excluded_ids ={}#



qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v9.db'#20mV run 10/09/25

all_run_ids= list(range(97, 138)) #v9 20mV on gate at instr; g4 dot;runs from 11 &13 /09/25
excluded_ids ={}#

all_run_ids= list(range(138, 179)) #v9 5mV on gate at instr; g4 dot;runs from 11 &13 /09/25
excluded_ids ={}#

all_run_ids= list(range(179, 220)) #v9 5mV on gate at instr; g4 dot;runs from 11 &13 /09/25
excluded_ids ={}#

all_run_ids= list(range(1944, 1984)) #v9 100uV on gate2 at instr; g3 dot;runs from 11 &13 /09/25
excluded_ids ={}

all_run_ids= list(range(2870, 2951)) #v9 100uV on gate2 at instr; g3 dot;runs from 11 &13 /09/25
excluded_ids ={}

qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v11.db'#20mV run 10/09/25

all_run_ids= list(range(42, 67)) #v9 100uV on gate2 at instr; g3 dot;runs from 11 &13 /09/25


qc.config["core"]["db_location"]=".\Data\Raw_data\CD12_B5_F4v41_01_12_25.db"#5mV 01/12/25

all_run_ids= list(range(6, 29))

all_run_ids= list(range(37, 90))


excluded_ids ={}

run_ids = [rid for rid in all_run_ids if rid not in excluded_ids]



freq,test_sweep=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)
sweeps = np.zeros((len(test_sweep), len(run_ids)))

sitpos_list=[]

for i, run_id in enumerate(run_ids):
   exp_name,_, sweep = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
   print(f"extracted_runid {run_id}")
   sweeps[:, i] = sweep
   sitpos=get_metadata(run_id,print_it=False,return_data=True)['qdac_ch06_dc_constant_V']
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
#plt.title("sweeps vs sitpos,-sb,312uV@instrg2(v9;runs 1261-1292)")
#plt.title("sweeps vs sitpos,-sb,20mV@instrg2g4dot(v9;runs 97-137)")
#plt.title("sweeps vs sitpos,-sb,5mV@instrg2g4dot(v9;runs 138-178)")
#plt.title("sweeps vs sitpos,-sb,1.25mV@instrg2g4dot(v9;runs 179-219)")
plt.title("sweeps vs sitpos,-sb,5mV@instr,28 holes(v41;runs 6-xx)")
plt.title("sweeps vs sitpos,-sb,5mV@instr,24 holes(v41;runs 6-xx)")
plt.show()
