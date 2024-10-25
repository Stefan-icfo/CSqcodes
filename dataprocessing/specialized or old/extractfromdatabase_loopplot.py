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
logtraces=[]
Vgates=[]
ids=[1453,1466,1467,1468,1469]
barriers=[0.6,0.7,0.5,0.3,0]
for id in ids:
    dataset_temp=qc.load_by_id(id)
    df_temp=dataset_temp.to_pandas_dataframe_dict()
#plot

    interdeps = dataset_temp.description.interdeps
    param_spec = interdeps.non_dependencies[0]  # hall resistance data
    #param_name = param_spec.name
    data_x = dataset_temp.get_parameter_data(param_spec)


    V_gate = np.array(data_x["G"]['QDAC_ch06_dc_constant_V'])

    trace=np.array(df_temp["G"])
    logtrace=np.log(trace)
    traces.append(trace)
    logtraces.append(logtrace)
    Vgates.append(V_gate)

#dataset_temp=qc.load_by_id(1871)
#df_temp=dataset_temp.to_pandas_dataframe_dict()
#trace2=np.array(df_temp["V_r"])

#dataset_temp=qc.load_by_id(1872)
#df_temp=dataset_temp.to_pandas_dataframe_dict()
#trace3=np.array(df_temp["V_r"])
#trace=trace3
#currentforward=np.array(df_temp["G"])
#currentbackward=np.array(df_1791["Current"])
#currentbackward=np.flip(currentbackward)
#smoothforward=centered_moving_average(trace,30)
#smoothforward=moving_average(trace,40)
#smoothbackward=np.flip(smoothforward)
#smooth=smoothforward

#frequency_forward=np.linspace(20e6,110e6,len(smoothforward))
#frequency_backward=np.linspace(10e6,200e6,len(smoothbackward))
for logtrace,Vgate,barrier in zip(logtraces,Vgates,barriers):
    plt.plot(Vgate,logtrace*1e6,label=f"{barrier} V")

#plt.plot(frequency_backward,smoothbackward)
#plt.plot(Vgates[0],traces[0],label='0.6')
#plt.plot(Vgate,traces[1],label='0.7')
    plt.xlabel("V_cs [V]")
    plt.ylabel("Conductance [uS]")
#plt.title("CD08_chipG7_devF6 5gmechanic 30smooth meas1871")
#plt.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
    plt.legend(loc="upper right")
    plt.show()

#example: make list with database entry
#voltagelist=df_1310["voltage"]