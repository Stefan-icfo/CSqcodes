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

dataset_temp15=qc.load_by_id(1240)
df_temp=dataset_temp15.to_pandas_dataframe_dict()
#plot
trace15=np.array(df_temp["I_rf"])
smoothforward15=moving_average(trace15,100)

dataset_temp5=qc.load_by_id(1089)
df_temp=dataset_temp5.to_pandas_dataframe_dict()
trace5=np.array(df_temp["Conductance"])
smoothforward5=moving_average(trace5,60)
#dataset_temp3=qc.load_by_id(891)
#df_temp=dataset_temp3.to_pandas_dataframe_dict()
#trace3=np.array(df_temp["I_rf"])
#smoothforward3=moving_average(trace3,120)
#dataset_temp10=qc.load_by_id(896)
#df_temp=dataset_temp10.to_pandas_dataframe_dict()
#trace10=np.array(df_temp["I_rf"])
#smoothforward10=moving_average(trace10,120)
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

frequency_forward=np.linspace(250e6,260e6,len(smoothforward15))

frequency_forward2=np.linspace(200e6,150e6,len(smoothforward5))

plt.plot(frequency_forward,smoothforward15,'b')
#plt.plot(frequency_forward2,smoothforward5,'r')
#plt.plot(frequency_forward,smoothforward5,'go')
#plt.plot(frequency_forward,smoothforward3,'ro')
#plt.plot(frequency_forward,smoothforward10,'yo')
plt.xlabel("Frequency(Hz)")
plt.ylabel("i")
#plt.title("mechan")
#plt.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
plt.show()

#example: make list with database entry
#voltagelist=df_1310["voltage"]