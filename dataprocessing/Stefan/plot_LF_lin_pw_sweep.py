import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re







qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v8.db'#20mV run 10/09/25

data_name='I_rf'
setpoint_name="zurich_oscs0_freq"

run_ids_list=[]
run_ids= list(range(1961, 1970)) #20u
run_ids_list.append(run_ids)
run_ids= list(range(1951, 1960)) #40u
run_ids_list.append(run_ids)
run_ids= list(range(1941, 1950)) #60u
run_ids_list.append(run_ids)
run_ids= list(range(1931, 1940)) #80u
run_ids_list.append(run_ids)
run_ids= list(range(1921, 1930)) #100u
run_ids_list.append(run_ids)
run_ids= list(range(1911, 1920)) #120u
run_ids_list.append(run_ids)
run_ids= list(range(1901, 1910)) #140u
run_ids_list.append(run_ids)
run_ids= list(range(1891, 1900)) #160u
run_ids_list.append(run_ids)
run_ids= list(range(1881, 1890)) #180u
run_ids_list.append(run_ids)
run_ids=  [1862,1863,1864,1865,1866,1873,1874,1875,1876,1877]#200uu
run_ids_list.append(run_ids)

avg_num=101


#freq,test_sweep=extract_1d(run_ids[0], data_1d_name = data_name, setpoint_name = setpoint_name,  plot = True)
#sweeps = np.zeros((len(test_sweep), len(run_ids)))

power_list=[]
maxf_list=[]

sweeps=[]
labels=["20uV","40uV","60uV","80uV","100uV","120uV","140uV","160uV","180uV","200uV"]
for runs,label_ in zip(run_ids_list,labels):
    for i, run_id in enumerate(runs):
        exp_name,freq, sweep = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
        print(f"extracted_runid {run_id}")
        sweeps.append(centered_moving_average(sweep,avg_num))
        sweepsum = np.sum(sweeps, axis=0)

    plt.plot(freq,sweepsum,label=label_)




#qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v9.db'#0mV run 
#runs= list(range(4, 13))
#for i, run_id in enumerate(runs):
#        exp_name,freq, sweep = extract_1d(run_id, data_1d_name=data_name, setpoint_name=setpoint_name, plot=False,return_exp_name=True)
#        print(f"extracted_runid {run_id}")
#        sweeps.append(centered_moving_average(sweep,avg_num))
#        sweepsum = np.sum(sweeps, axis=0)

#plt.plot(freq,sweepsum,label="0uV")
plt.xlabel("Frequency [MHz]")
plt.ylabel("current")
plt.legend()


plt.title("powersweep_g2dot")

plt.show()

