import matplotlib.pyplot as plt
from dataprocessing.extract_fkts_save import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re







data_name_1='G'
data_name_2='I_sens'
setpoint_name="QDAC_ch06_dc_constant_V"



qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v32_19_11_25.db'




all_run_ids= list(range(13, 54)) #Vg_cs_up
excluded_ids ={38}
run_ids = [rid for rid in all_run_ids if rid not in excluded_ids]
print(f"{len(run_ids)} run_ids")

#all_run_ids= list(range(403, 443)) ##Vg_cs_down
#excluded_ids ={}
#run_ids = [rid for rid in all_run_ids if rid not in excluded_ids]
#print(f"{len(run_ids)} run_ids")




Vg_cs,G_test=extract_1d(run_ids[0], data_1d_name = data_name_1, setpoint_name = setpoint_name,  plot = True)
#Bias = np.arange(-100, 301, 10)

Bias = np.arange(-100, 291, 10)
print(f"{len(Bias)} Bias values")
map_G = np.zeros((len(G_test), len(run_ids)))
map_Isens = np.zeros((len(G_test), len(run_ids)))


for i, run_id in enumerate(run_ids):
   #fill initial sweep
   exp_name,_, G = extract_1d(run_id, data_1d_name=data_name_1, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
   exp_name,_, I_sens = extract_1d(run_id, data_1d_name=data_name_2, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
   print(f"extracted_runid {run_id} of first run at this Voltage")
   map_G[:, i] = G 
   map_Isens[:, i] = I_sens

    

#plt.pcolormesh(freq[::-1] / 1e6, drives, spectra.T, cmap='magma', shading='auto')  # reverse frequency axis for other sideband
mesh = plt.pcolormesh(Vg_cs,Bias,map_G.T,cmap='magma',shading='auto')

# Set color limits on the returned QuadMesh object
#mesh.set_clim(0,2)

plt.xlabel("Vg_cs [V]")
plt.ylabel("Bias [uV]")
plt.colorbar(label="G [S]")
plt.show()

mesh = plt.pcolormesh(Vg_cs,Bias,map_Isens.T,cmap='magma',shading='auto')

# Set color limits on the returned QuadMesh object
#mesh.set_clim(0,2)

plt.xlabel("Vg_cs [V]")
plt.ylabel("Bias [uV]")
plt.colorbar(label="I_sens [A]")
plt.show()