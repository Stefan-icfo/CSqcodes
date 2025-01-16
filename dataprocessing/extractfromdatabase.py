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
qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1_part2.db'

experiments=qc.experiments()

##dataset_temp=qc.load_by_id(1789)
#df_1789=dataset_temp.to_pandas_dataframe_dict()

dataset_temp=qc.load_by_id(1134)
df_temp=dataset_temp.to_pandas_dataframe_dict()
#plot

interdeps = dataset_temp.description.interdeps
param_spec = interdeps.non_dependencies[0]  # hall resistance data
#param_name = param_spec.name
data_x = dataset_temp.get_parameter_data(param_spec)


    #V_gate = np.array(data_x["Conductance"]['QDAC_ch01_dc_constant_V'])

trace=np.array(df_temp["I_rf"])
num_points = len(trace)  # Get the number of points in v_r data
time_array = np.linspace(2, 5, num_points)  # Create time array from 2 to 5 seconds



plt.plot(time_array,trace)

plt.show()

