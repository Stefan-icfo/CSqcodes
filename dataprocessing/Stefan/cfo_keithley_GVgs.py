import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re







qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v9.db'#20mV run 10/09/25

data_name='Conductance_IV'
setpoint_name="QDAC_ch02_dc_constant_V"

run_ids_nodrive=[30,31,32,33,34,36,38,40,42]
run_ids_drive750=[35,37,39,41,43,47,48,49,50]



run_ids_drive375=list(range(51,60))
#labels=["25","27"]


#freq,test_sweep=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)
#sweeps = np.zeros((len(test_sweep), len(run_ids)))

nodrive_G_IV_list=[]

for run_id in run_ids_nodrive:
        exp_name,G2V, G_IV = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
        print(f"extracted_runid {run_id}")
        nodrive_G_IV_list.append(G_IV)

nodrive_sum = np.sum(nodrive_G_IV_list, axis=0)    
plt.plot(G2V,nodrive_sum,label="nodrive")

integral_nodrive=sum(nodrive_sum)


drive_G_IV_list750=[]
for run_id in run_ids_drive750:
        exp_name,G2V, G_IV = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
        print(f"extracted_runid {run_id}")
        drive_G_IV_list750.append(G_IV)

drive_sum750 = np.sum(drive_G_IV_list750, axis=0)    
plt.plot(G2V-250e-6,drive_sum750*672/538,label="drive 750mV")
scaled_integral_750m=sum(drive_sum750*672/538)

drive_G_IV_list375=[]
for run_id in run_ids_drive375:
        exp_name,G2V, G_IV = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
        print(f"extracted_runid {run_id}")
        drive_G_IV_list375.append(G_IV)

drive_sum375= np.sum(drive_G_IV_list375, axis=0)    
plt.plot(G2V-200e-6,drive_sum375*672/617,label="drive 375mV")
scaled_integral_375m=sum(drive_sum375*672/617)




#qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v9.db'#0mV run 
#runs= list(range(4, 13))
#for i, run_id in enumerate(runs):
#        exp_name,freq, sweep = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
#        print(f"extracted_runid {run_id}")
#        sweeps.append(centered_moving_average(sweep,avg_num))
#        sweepsum = np.sum(sweeps, axis=0)

#plt.plot(freq,sweepsum,label="0uV")
plt.xlabel("G2_Voltage")
plt.ylabel("IV_conductance")
plt.legend()


plt.title("calibrate_G2_rf")

plt.show()


plt.plot(G2V,nodrive_sum,label="nodrive")
plt.plot(G2V,drive_sum750,label="drive 750mV")
plt.plot(G2V,drive_sum375,label="drive 375mV")
plt.xlabel("G2_Voltage")
plt.ylabel("IV_conductance")
plt.legend()


plt.title("calibrate_G2_rf")

plt.show()

plt.plot([0,375,750],[integral_nodrive,scaled_integral_375m,scaled_integral_750m],'*')
#plt.ylim(0,0.4e-3)
x_vals = np.array([0, 375, 750])
y_vals = np.array([integral_nodrive, scaled_integral_375m, scaled_integral_750m])
plt.plot(x_vals, np.poly1d(np.polyfit(x_vals, y_vals, 1))(x_vals), '-r')
plt.show()