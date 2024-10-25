import math
import numpy as np
#import math
#import csv
import matplotlib.pyplot as plt
#import os


import pandas as pd
import qcodes as qc

def moving_average(a, n=30):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

#database location
qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'QuantumSimulator.db'

experiments=qc.experiments()

##dataset_temp=qc.load_by_id(1789)
#df_1789=dataset_temp.to_pandas_dataframe_dict()
for n in range(9): #1798+16 #1814+16 #1830+9
    m=1878+n
    dataset_temp=qc.load_by_id(m)
    df_temp=dataset_temp.to_pandas_dataframe_dict()
    smoothtrace=np.array(df_temp["V_r"])
    #currentbackward=np.array(df_1791["Current"])
    #currentbackward=np.flip(currentbackward)
    smoothtrace=moving_average(smoothtrace,30)
    #smoothbackward=moving_average(currentbackward,10)


    frequency=np.linspace(75e6,100e6,len(smoothtrace))
    plt.plot(frequency,smoothtrace)
#df_1140=dataset_1140.to_pandas_dataframe_dict()
#df_1141=dataset_1141.to_pandas_datafrae_dict()

plt.show()
#plot

#currentforward=np.array(df_1838["Current"])
#currentbackward=np.array(df_1791["Current"])
#currentbackward=np.flip(currentbackward)
#smoothforward=moving_average(currentforward,10)
#smoothbackward=moving_average(currentbackward,10)


#frequency_forward=np.linspace(80e6,100e6,len(smoothforward))
#frequency_backward=np.linspace(10e6,200e6,len(smoothbackward))

#plt.plot(frequency_forward,smoothforward)
#plt.plot(frequency_backward,smoothbackward)

#plt.xlabel("voltage")
#plt.ylabel("current")
#plt.title("CD08_chipG7_devF6 4gmechanicspowe-2volt779")
#plt.grid(which='major', linestyle=':', linewidth='0.5', color='gray')


#example: make list with database entry
#voltagelist=df_1310["voltage"]