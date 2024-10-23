import math
import numpy as np
#import math
#import csv
import matplotlib.pyplot as plt
#import os


import pandas as pd
import qcodes as qc
from utils.CS_utils import centered_moving_average

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

#database location
qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1.db'

experiments=qc.experiments()

##dataset_temp=qc.load_by_id(1789)
#df_1789=dataset_temp.to_pandas_dataframe_dict()
traces=[]
ids=[1391,1392,1393,1394,1398,1403,1404]
temps=[48,41,36,33,30,28,23]
for id in ids:
    dataset_temp=qc.load_by_id(id)
    df_temp=dataset_temp.to_pandas_dataframe_dict()
#plot

    interdeps = dataset_temp.description.interdeps
    param_spec = interdeps.non_dependencies[0]  # hall resistance data
    #param_name = param_spec.name
    data_x = dataset_temp.get_parameter_data(param_spec)


    V_gate = np.array(data_x["Conductance"]['QDAC_ch01_dc_constant_V'])

    trace=np.array(df_temp["Conductance"])
    traces.append(trace)



#frequency_forward=np.linspace(20e6,110e6,len(smoothforward))
#frequency_backward=np.linspace(10e6,200e6,len(smoothbackward))
for trace,temp in zip(traces,temps):
    plt.plot(V_gate,trace*1e6,label=f"{temp} K")

#plt.plot(frequency_backward,smoothbackward)

plt.xlabel("V_cs [V]")
plt.ylabel("Conductance [uS]")
#plt.title("CD08_chipG7_devF6 5gmechanic 30smooth meas1871")
#plt.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
plt.legend(loc="upper right")
plt.show()

#example: make list with database entry
#voltagelist=df_1310["voltage"]