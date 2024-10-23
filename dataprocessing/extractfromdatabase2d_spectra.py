import math
#import math
#import csv
import matplotlib.pyplot as plt
#import os


import pandas as pd
import qcodes as qc
import numpy as np

def voltage_to_dbm(v_rms, rbw, impedance=50):
    # Convert voltage to dBm using the formula
    dbm = 20 * np.log10(v_rms) + 10 * np.log10(1 / (impedance * rbw)) + 30
    return dbm
def voltage_to_psd(v_rms, rbw, impedance=50):
  
    # Calculate PSD using the formula
    psd = (v_rms ** 2) / (impedance * rbw)
    return psd

#database location
qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1_zurichdata.db'
experiments=qc.experiments()

dataset=qc.load_by_id(97)


pdf=dataset.to_pandas_dataframe_dict()
freq_spec_raw=pdf["Voltage_fft_avg"]
freq_spec_np=np.array(freq_spec_raw)
# ---------------------Geting the data from the database---------------------
# pprint(dataset.get_parameter_data())
interdeps = dataset.description.interdeps
param_spec = interdeps.non_dependencies[0]  # hall resistance data
param_name = param_spec.name
data_xy = dataset.get_parameter_data(param_spec)
xy = data_xy[param_name][param_name]

#g1:outer gate
#g2:inner gate

time_raw = data_xy[param_name]['time_param']
freq_raw = data_xy[param_name]['freq_param']
time_np=np.array(time_raw)
freq_np=np.array(freq_raw)
time=np.unique(time_np)
freq=np.unique(freq_np)
#real_time_len
nr_time_points=round(len(freq_spec_np)/len(freq))
real_time_axis=np.linspace(0,4.772*nr_time_points,nr_time_points)

freq_spectrum_real=np.zeros([nr_time_points, len(freq)])


for m in range(nr_time_points):
    for n in range(len(freq)):
        freq_spectrum_real[m,n]=freq_spec_np[m*len(freq)+n]

        #integral


plt.pcolor(freq,real_time_axis,freq_spectrum_real)
plt.xlabel("frequency delta from demod [Hz]")
plt.ylabel("time [s]")
plt.show()

average=True
if average:
    time_avg=np.mean(freq_spectrum_real,axis=0)
    plt.plot(freq,time_avg)
    plt.xlabel("frequency delta from demod [Hz]")
    plt.ylabel("Voltage [V]")
    plt.show()

    plt.plot(freq,voltage_to_psd(time_avg,0.209))
    plt.xlabel("frequency delta from demod [Hz]")
    plt.ylabel("Power [W/Hz]")
    #plt.ylim(top=5.5e-13)
    plt.show()


