import math
import numpy as np
#import math
#import csv
import matplotlib.pyplot as plt
#import os


import pandas as pd
import qcodes as qc

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

#database location
qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'QuantumSimulator.db'

experiments=qc.experiments()

##dataset_temp=qc.load_by_id(1789)
#df_1789=dataset_temp.to_pandas_dataframe_dict()

dataset_temp=qc.load_by_id(1878)
df_temp=dataset_temp.to_pandas_dataframe_dict()
#plot
trace=np.array(df_temp["V_r"])

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
smoothforward=moving_average(trace,30)
#smoothbackward=moving_average(currentbackward,10)
smoothbackward=np.flip(smoothforward)
#smooth=smoothforward

frequency_forward=np.linspace(82e6,86e6,len(smoothbackward))
#frequency_backward=np.linspace(10e6,200e6,len(smoothbackward))

#Lorentzian
Gamma1=36e3
Fresonance1=84.609e6
Amax1=3.85e-6 #peak of Lortenzian, in dBm
Background1=3.7e-6 #in a.u.


Fnp=frequency_forward[2500:2680]
LorentzianLin=Amax1*(((Gamma1/2)**2)/((Fresonance1-Fnp)**2+(Gamma1/2)**2))+Background1



plt.plot(frequency_forward[2500:2680],smoothbackward[2500:2680])
plt.plot(frequency_forward[2500:2680],LorentzianLin, linewidth=3)
#plt.plot(frequency_backward,smoothbackward)



LorentzianLin=(((Gamma1/2)**2)/((Fresonance1-Fnp)**2+(Gamma1/2)**2))+Background1

plt.xlabel("Frequecy(Hz)")
plt.ylabel("Amplitude (a.u.)")
#plt.title("CD08_chipG7_devF6 5gmechanic 30smooth meas1871")
#plt.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
plt.show()

#example: make list with database entry
#voltagelist=df_1310["voltage"]